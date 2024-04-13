% runAssignment5_Group16
%  group 16, AY2023-2024
% Compute the 
%
% to run:
% > runAssignment5_Group16

% clear workspace
clear all;
close all;
clc;

%% Settings

formatData ='dd/MM/yyyy'; % Pay attention to your computer settings 

rng(42);   % Fix the random number generator ("the answer to Life, the Universe, and Everything")

tic;       % Set the clock to find the time of execution

%% Data import

% import the discounts factors
load("discounts.mat");
% load the data for the second and third exercises
load("cSelect20230131_B.mat");

%% Point 1: Certificate Pricing

% spread over libor
s_spol = 100 * 1e-4; % 100 bps
% upfront percentage of notional
X = 2 / 100; % 2%
P = 95 / 100; % 95%
principal_Amount = 100e6; % 100M
% no counterparty risk and neglect IR dynamics

% parameters for the underlying (Black model)

% initial date (datetime) 18/02/2008
t0 = datetime(2008, 2, 18, 'Format', formatData);
% initial prices
ENEL_0 = 100;
AXA_0 = 200;
% volatilities
sigma_ENEL = 16.1 / 100; % 16.1%
sigma_AXA = 20 / 100; % 20%
rho = 40 / 100; % 40%
% dividend yields
d_ENEL = 2.5 / 100; % 2.5%
d_AXA = 2.7 / 100; % 2.7%

%% Find the dates of the cash flows

% party B dates (annualy for 4 years), 4 dates
partyB_dates = t0 + calyears(1:4);
% move to business days if needed (modified follow convention, no holidays)
partyB_dates(~isbusday(partyB_dates,0)) = busdate(partyB_dates(~isbusday(partyB_dates,0)), "modifiedfollow", 0);
disp('Party B dates: ');
disp(partyB_dates);
% convert to datenums
partyB_dates = datenum(partyB_dates);

% party A dates (quarterly for 4 years) 12 dates
partyA_dates = t0 + calmonths(3:3:48);
partyA_dates(~isbusday(partyA_dates,0)) = busdate(partyA_dates(~isbusday(partyA_dates,0)), "follow", 0);
disp('Party A dates: ');
disp(partyA_dates);
partyA_dates = datenum(partyA_dates);

%% Compute the partecipation coefficient

N_sim = 1e7;

alpha = priceCertificate(ENEL_0,  sigma_ENEL, d_ENEL, AXA_0, sigma_AXA, d_AXA, rho, s_spol, P, X, ...
    principal_Amount, N_sim, partyA_dates, partyB_dates, dates, discounts);

disp(['The participation coefficient is: ', num2str(alpha)]);

%% Point 2: Pricing Digital Option

% Price with Black Formula
Notional = 1e7;
% strike
k = cSelect.reference;
% maturity
ACT_365 = 3;
T = yearfrac(t0, t0 + calyears(1), ACT_365);
% value of the digital option at maturity if s > k
payment = 0.05 * Notional;
% compute the discount factor at 1 year
discount_1y = intExtDF(discounts, dates, datenum(t0 + calyears(1)));
% find the corrisponding interest rate
r = -log(discount_1y) / T;
S_digital  = k * discount_1y;

% Load volatility smile
strikes = cSelect.strikes;
surface = cSelect.surface;

sigma_digial = interp1(strikes, surface, k, 'spline');

plot(strikes, surface); hold on;
plot(k, sigma_digial,'x')

d_1 = (log(k / k) + (0.5 * sigma_digial^2) * T) / (sigma_digial * sqrt(T));
d_2 = d_1 - sigma_digial * sqrt(T);

price_digital_black = payment * discount_1y * normcdf(d_2)

% Now use a volatility approach
% Find k_1 and k_2 that contain k
k_1 = strikes(find(strikes < k, 1, 'last'));
k_2 = strikes(find(strikes > k, 1, 'first'));

% Find the corresponding volatilities
sigma_1 = interp1(strikes, surface, k_1, 'spline');
sigma_2 = interp1(strikes, surface, k_2, 'spline');

% compute the skew in that point
m = (sigma_2 - sigma_1) / (k_2 - k_1);

% Compute the vega under black model
vega = S_digital * normpdf(d_1) * sqrt(T);

% Now compute the digital price

price_digital_implied = price_digital_black - vega * m

% now implement monte carlo simulation
N = 1e7;
Z = randn(N, 1);
S = S_digital * exp((r - 0.5 * sigma_digial^2) * T + sigma_digial * sqrt(T) * Z);
% payoff digital option pays 0.05 
payoff = payment *  (S > k);
price_digital_monte_carlo = mean(payoff) * discount_1y

% compute the error
error = abs(price_digital_implied - price_digital_black);
disp(['The error between the implied and black price is: ', num2str(error*1e4), ' bps']);

%% Point 3: Pricing

% parameters of the mean-variance mixture model
alpha = 0.5;
sigma = 20 / 100;
kappa = 1;
eta = 3;
t = 1;
% moneyness
x = (-25:1:25) / 100;
S_0 = cSelect.reference;
F_0 = S_0 / discount_1y;

%% Point 3.a: FFT method, alpha = 1/2

% compute the call prices with the FFT method
M_FFT = 15;
flag = 'FFT';
callPrices_FFT = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, flag);

% put the results in a txt file
fileID = fopen('callPrices_FFT.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_FFT);
fclose(fileID);

%% Point 3.b: Quadrature method, alpha = 1/2
% compute the call prices with the quadrature method
M_quad = 21;
flag = 'quad';
callPrices_quad = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_quad, flag);

% put the results in a txt file
fileID = fopen('callPrices_quad.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_quad);
fclose(fileID);

%% Point 3.c: Monte Carlo simulation with alpha = 1/2

% use a MonteCarlo simulation to compute the call prices
N = 1e7;

% compute the Laplace exponent
ln_L = @(omega) t/kappa * (1 - alpha)/alpha * ...
    (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );
% draw the standard normal random variables
g = randn(N, 1);
% draw the inverse gaussian random variables
G = random('inversegaussian', 1, kappa/t, N, 1);

ft = sqrt(t) * sigma * sqrt(G) .* g - (0.5 + eta) * t * sigma^2 * G - ln_L(eta);

FT = F_0 * exp(ft);

% compute the call prices
callPrices_MC = zeros(size(x));
for i = 1:length(x)
    callPrices_MC(i) = mean(max(FT - exp(-x(i))*F_0, 0)) * discount_1y;
end

% put the results in a txt file
fileID = fopen('callPrices_MC.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_MC);
fclose(fileID);

%% Point 3.c: Black prices (check)

% compute the real price
realVols = cSelect.surface;
realStrikes = cSelect.strikes;
realDividend = cSelect.dividend;
% for each strike compute the black price
realPrices = zeros(size(realStrikes));
for i = 1:length(realStrikes)
    d_1 = (log(F_0 / realStrikes(i)) + (-realDividend + 0.5 * realVols(i)^2) * t) / (realVols(i) * sqrt(t));
    d_2 = d_1 - realVols(i) * sqrt(t);
    realPrices(i) = discount_1y * F_0 * exp(-realDividend * t) * normcdf(d_1) - realStrikes(i) * discount_1y * normcdf(d_2);
end

%% Plot results with alpha = 1/2

figure;

% plot quadrature
plot(x, callPrices_quad);
hold on
% plot FFT
plot(x, callPrices_FFT);
% plot the Monte Carlo
plot(x, callPrices_MC);
% plot the Black prices
hold on
plot(log(F_0 ./ realStrikes), realPrices, 'x');

title('Call prices with different methods and alpha = 1/2');
xlabel('Moneyness');
legend('Quadrature', 'FFT', 'Monte Carlo', 'Black prices');

%% Point 3.d: Use alpha = 2/3

% run the FFT with alpha= 2/3
alpha = 2/3;
callPrices_FFT_2_3 = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, 'FFT');
callPrices_quad_2_3 = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, 'quad');

% put the results in a txt file
fileID = fopen('callPrices_FFT_2_3.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_FFT_2_3);
fclose(fileID);

%% Point 3.d: Plot the results

% plot the call prices
figure;
% plot FFT
plot(x, callPrices_FFT_2_3);
hold on
% plot quadrature
plot(x, callPrices_quad_2_3);
% plot old results
plot(x, callPrices_FFT);
title('Call prices with different methods and alpha = 2/3');
xlabel('Moneyness');
legend('FFT', 'Quadrature', 'FFT alpha = 1/2');

%% Point 4: Volatility Surface Calibration

% alpha = 1/3
alpha = 1/3;
% compute the log moneyess from the strikes
log_moneyness = log(realStrikes / F_0);

% create a function that the prices of the call options given the strikes
prices = @(sigma, kappa, eta) callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, log_moneyness, M_FFT, 'FFT');

% compute the lower bound for eta
omega_down = (1 - alpha) / (kappa * sigma^2);

% calibrate the model using fmincon
% initial guess
x0 = [0.2, 1, 1];

% lower bounds
lb = [0, 0, -omega_down];

% calibration
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4);

[x, fval] = fmincon(@(x) norm(prices(x(1), x(2), x(3)) - realPrices), x0, [], [], [], [], lb, [], [], options);

% compute the prices with the calibrated parameters
prices_calibrated = prices(x(1), x(2), x(3));

% invert the prices using black formula
implied_vols = zeros(size(realStrikes));
for i = 1:length(realStrikes)
    d_1 = (log(F_0 / realStrikes(i)) + (-realDividend + 0.5 * implied_vols(i)^2) * t) / (implied_vols(i) * sqrt(t));
    d_2 = d_1 - implied_vols(i) * sqrt(t);
    implied_vols(i) = fmincon(@(x) abs(prices(x(1), x(2), x(3)) - realPrices(i)), x0, [], [], [], [], lb, [], [], options);
end

% plot the results
figure;
plot(realStrikes, implied_vols);
hold on
plot(realStrikes, realVols, 'x');
title('Implied volatilities');
xlabel('Strikes');
legend('Implied volatilities', 'Real volatilities');

toc