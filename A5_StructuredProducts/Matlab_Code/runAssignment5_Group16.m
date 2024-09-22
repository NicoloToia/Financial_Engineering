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
close all;

%% Settings

addpath Data\
addpath Functions\

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
t0 = datetime(dates(1), 'ConvertFrom', 'datenum');
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
% convert to datenums
partyB_dates = datenum(partyB_dates);

% party A dates (quarterly for 4 years) 12 dates
partyA_dates = t0 + calmonths(3:3:48);
partyA_dates(~isbusday(partyA_dates,0)) = busdate(partyA_dates(~isbusday(partyA_dates,0)), "follow", 0);
partyA_dates = datenum(partyA_dates);

%% Compute the partecipation coefficient

% set number of simulations
N_sim = 1e7;

% compute the participation coefficient
[alpha, IC_alpha] = priceCertificate(ENEL_0,  sigma_ENEL, d_ENEL, AXA_0, sigma_AXA, d_AXA, rho, s_spol, P, X, ...
    principal_Amount, N_sim, partyA_dates, partyB_dates, dates, discounts, 0.95);

% display the results
disp(['The participation coefficient is: ', num2str(alpha)]);
disp(['The confidence interval is: [', num2str(IC_alpha(1)), ', ', num2str(IC_alpha(2)), ']']);

% put the results in a txt file
fileID = fopen('patrtecipation_coefficient.txt', 'w');
fprintf(fileID, '%12.8f\n', alpha);
fprintf(fileID, '%12.8f\n', IC_alpha);
fclose(fileID);


%% Point 2: Pricing Digital Option

% Price with Black Formula
Notional = 1e7;
% spot and strike from data
S_0 = cSelect.reference;
k = S_0;
d = cSelect.dividend;
% maturity
ACT_365 = 3;
T = yearfrac(t0, t0 + calyears(1), ACT_365);
% compute the discount factor at 1 year
discount_1y = intExtDF(discounts, dates, datenum(t0 + calyears(1)));
% find the corrisponding zero rate
r = -log(discount_1y) / T;
% compute the forward price
F_0  = S_0 / discount_1y * exp(-d * T);

% Load volatility smile
strikes = cSelect.strikes;
surface = cSelect.surface;
% sigma digital
sigma_digital = interp1(strikes, surface, k, 'spline');
% plot the volatility smile
plot(strikes, surface);
hold on;
plot(k, sigma_digital,'x', 'MarkerSize', 5, 'LineWidth', 5);
legend('Volatility smile', 'Volatility at the money');

% compute the price via DigitalPrice function
% flag = 1: Black formula
% flag = 2: Volatility Approach
% flag = 3: Monte Carlo (Black dynamics)
price_digital_black = DigitalPrice(Notional , T , F_0 , d , discount_1y , k , strikes , surface , 1);
price_digital_implied = DigitalPrice(Notional , T , F_0 , d , discount_1y , k , strikes , surface , 2);
price_digital_monte_carlo = DigitalPrice(Notional , T , F_0 , d , discount_1y , k , strikes , surface , 3);

% compute the error
error = abs(price_digital_implied - price_digital_black);

error_percentage = (error / price_digital_black) * 100;

% display results
fprintf(['\nThe black price is: ', num2str(price_digital_black)]);
fprintf(['\nThe implied price is: ', num2str(price_digital_implied)]);
fprintf('\nThe error between the implied and black price is: %.2f€ which is %.2f%% of the black price', error, error_percentage);
fprintf(['\nThe monte carlo price is: ', num2str(price_digital_monte_carlo), '\n']);

fileID = fopen('digital_prices.txt', 'w');
fprintf(fileID, '%12.8f\n', price_digital_black);
fprintf(fileID, '%12.8f\n', price_digital_implied);
fprintf(fileID, '%12.8f\n', price_digital_monte_carlo);
fclose(fileID);

plot_payoff = plotpayoff(strikes, k);

%% Point 3: Pricing

% parameters of the mean-variance mixture model
alpha = 0.5;
sigma = 20 / 100;
kappa = 1;
eta = 3;
t = 1;
% moneyness
x = (-25:1:25) ./ 100;
S_0 = cSelect.reference;
d = cSelect.dividend;
F_0 = S_0 / discount_1y * exp(-d * t);

%% Point 3.a: FFT method, alpha = 1/2

% compute the call prices with the FFT method
M_FFT = 15;
dz = 0.0025;
flag = 'FFT';
callPrices_FFT = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, dz, flag);

% put the results in a txt file
fileID = fopen('callPrices_FFT.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_FFT);
fclose(fileID);

%% Point 3.b: Quadrature method, alpha = 1/2
% compute the call prices with the quadrature method
M_quad = 20;
flag = 'quad';
callPrices_quad = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_quad, dz, flag);

% put the results in a txt file
fileID = fopen('callPrices_quad.txt', 'w');
fprintf(fileID, '%12.8f\n', callPrices_quad);
fclose(fileID);

% compute the MSE between the two methods
distance = 1 / length(x) * sum((callPrices_FFT - callPrices_quad).^2);

if distance < 1e-4 % one bp of tolerance
    disp('The two methods converge');
else
    disp('The two methods do not converge');
end

%% Point 3.c: Monte Carlo simulation with alpha = 1/2

% use a MonteCarlo simulation to compute the call prices
N = 1e7;

% compute the Laplace exponent
ln_L = @(omega) t/kappa * (1 - alpha)/alpha * ...
    (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );
% draw the standard normal random variables
g = randn(N, 1);
% draw the inverse gaussian random variables
G = random('inversegaussian', 1, t/kappa, N, 1);

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

% check moments of inverse gaussian

% numerical first 4 moments
mu_1 = mean(G);
mu_2 = mean(G.^2);
mu_3 = mean(G.^3);
mu_4 = mean(G.^4);

% analytical first 4 moments
mu = 1;
lambda = t/kappa;
mu_1_an = 1;
mu_2_an = mu^2 * (lambda + mu) / lambda;
mu_3_an = mu^3 * (lambda^2 + 3 * lambda * mu + 3 * mu^2) / lambda^2;
mu_4_an = 5 * mu^2 / lambda * mu_3_an + mu_2_an * mu^2;

% print a table of the moments

disp('The first 4 moments of the inverse gaussian distribution are:');
disp(' ');
disp('Numerical | Analytical');
disp('-----------------------');
disp([mu_1, mu_1_an]);
disp([mu_2, mu_2_an]);
disp([mu_3, mu_3_an]);
disp([mu_4, mu_4_an]);

%% Point 3.c: Black prices (check)

% compute the real price
realVols = cSelect.surface;
realStrikes = cSelect.strikes;
% for each strike compute the black price
realPrices = zeros(size(realStrikes));
for i = 1:length(realStrikes)
    realPrices(i) = blackPrice(F_0, realStrikes(i), t, realVols(i), discount_1y, 'call');
end

%% Plot results with alpha = 1/2

figure;

% plot quadrature
plot(x, callPrices_quad,'--x');
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
legend('Quadrature', 'FFT', 'Monte Carlo','Black prices');

%% Point 3.d: Use alpha = 2/3

% run the FFT with alpha= 2/3
alpha = 2/3;
M_FFT = 15;
callPrices_FFT_2_3 = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, dz, 'FFT');
callPrices_quad_2_3 = callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, x, M_FFT, dz, 'quad');

% compute the MSE between the two methods
distance = 1 / length(x) * sum((callPrices_FFT_2_3 - callPrices_quad_2_3).^2);

if distance < 1e-4 % one bp of tolerance
    disp('The two methods converge');
else
    disp('The two methods do not converge');
end

% compute the average difference between the two methods
average_diff = mean(abs(callPrices_FFT_2_3 - callPrices_FFT));
average_diff_pct = mean(abs(callPrices_FFT_2_3 - callPrices_FFT)) / mean(callPrices_FFT) * 100;
disp(['The average difference between the two alphas is: ', num2str(average_diff)]);
disp(['The average difference in percentage is: ', num2str(average_diff_pct)]);

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
log_moneyness = log(F_0 ./ realStrikes);

% create a function that the prices of the call options given the strikes
prices = @(sigma, kappa, eta) callIntegral(discount_1y, F_0, alpha, sigma, kappa, eta, t, log_moneyness, M_FFT, dz, 'FFT');

% compute the lower bound for eta
omega_down = (1 - alpha) / (kappa * sigma^2);

% calibrate the model using fmincon
% initial guess
x0 = [0.2, 1, 3];

% lower bounds
lb = [0, 0, -omega_down];

% calibration
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4, 'Display', 'off');

[x, fval] = fmincon(@(x) norm(prices(x(1), x(2), x(3)) - realPrices), x0, [], [], [], [], lb, [], [], options);

% display the results
disp(['Calibrated parameters']);
disp(['Sigma: ', num2str(x(1))]);
disp(['Kappa: ', num2str(x(2))]);
disp(['Eta: ', num2str(x(3))]);

% compute the prices with the calibrated parameters
prices_calibrated = prices(x(1), x(2), x(3));

% plot the results
figure;
plot(realStrikes, prices_calibrated);
hold on
plot(realStrikes, realPrices, 'x');
title('Calibrated prices');
xlabel('Strikes');
legend('Calibrated prices', 'Real prices');

%% Point 4: plot the model implied volatilities

% invert the prices using black formula
model_implied_vols = zeros(size(realStrikes));
for i = 1:length(realStrikes)
    callPrice = @(s) blackPrice(F_0, realStrikes(i), t, s, discount_1y, 'call');
    % initial guess and lower bound
    x0 = 0.2; lb = 0;
    % compute the price with the calibrated parameters
    model_price = prices_calibrated(i);
    model_implied_vols(i) = fmincon(@(s) abs(model_price - callPrice(s)), x0, [], [], [], [], lb, [], [], options);
end

% plot the results
figure;
plot(realStrikes, model_implied_vols);
hold on
plot(realStrikes, realVols, 'x');
title('Implied volatilities');
xlabel('Strikes');
legend('Implied volatilities', 'Real volatilities');

% compute the error between the two volatilities as MAPE
distance = 1 / length(realStrikes) * sum(abs(realVols - model_implied_vols) ./ realVols) * 100;

disp(['The MAPE between the two volatilities is: ', num2str(distance)]);

% compute the maximum percentage error
max_error = max(abs(realVols - model_implied_vols) ./ realVols) * 100;

disp(['The MRE is: ', num2str(max_error)]);

% compute computation time
disp(['The computation time is: ', num2str(toc), ' seconds']);