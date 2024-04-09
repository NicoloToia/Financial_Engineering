% runAssignment5_Group16
%  group 16, AY2023-2024
% Compute the 
%
% to run:
% > runAssignment2_Group16

% clear workspace
clear all;
close all;
clc;

% NOTE : for long computational time the number of iterations are set by
% default at 1e4, use a larger number if more precision is required

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
sigma_ENEL = 16.1 / 100;
sigma_AXA = 20 / 100;
rho = 40 / 100;
% dividend yields
d_ENEL = 2.5 / 100;
d_AXA = 2.7 / 100;

%% Find the dates of the cash flows

% fixed leg (annualy for 4 years)
fixedLegDates = t0 + calyears(0:4);
% move to business days if needed (follow convention, no holidays)
fixedLegDates(~isbusday(fixedLegDates,0)) = busdate(fixedLegDates(~isbusday(fixedLegDates,0)), "follow", 0);
% convert to datenums
fixedLegDates = datenum(fixedLegDates);

% floating leg (quarterly for 4 years)
floatingLegDates = t0 + calmonths(0:3:48);
floatingLegDates(~isbusday(floatingLegDates,0)) = busdate(floatingLegDates(~isbusday(floatingLegDates,0)), "follow", 0);
disp(floatingLegDates)
floatingLegDates = datenum(floatingLegDates);

%% Simulate underlying prices and compute the party B leg

% compute the year fraction between fixedLegDates
ACT_365 = 3;
deltas_fixed = yearfrac(fixedLegDates(1:end-1), fixedLegDates(2:end), ACT_365);

% compute the discount factors for the fixed leg
discounts_fixed = intExtDF(discounts, dates, fixedLegDates);
% compute the forward discount factors
forward_discounts_fixed = discounts_fixed(2:end) ./ discounts_fixed(1:end-1);

% number of simulations
N = 1e6;

% variables to store the simulations
ENEL = zeros(N, 5);
AXA = zeros(N, 5);
% set initial values
ENEL(:,1) = ENEL_0;
AXA(:,1) = AXA_0;

% for each time step
for i = 1:4
    % extract the two gaussian random variables
    Z = mvnrnd([0 0], [1 rho; rho 1], N);
    Z_ENEL = Z(:, 1);
    Z_AXA = Z(:, 2);
    % update the prices
    ENEL(:, i+1) = 1/forward_discounts_fixed(i) * ENEL(:, i) .* exp((-d_ENEL - 0.5 * sigma_ENEL^2) * deltas_fixed(i) + sigma_ENEL * sqrt(deltas_fixed(i)) * Z_ENEL);
    AXA(:, i+1) = 1/forward_discounts_fixed(i) * AXA(:, i) .* exp((-d_AXA - 0.5 * sigma_AXA^2) * deltas_fixed(i) + sigma_AXA * sqrt(deltas_fixed(i)) * Z_AXA);
end 

% compute the returns each year
ratio_ENEL = ENEL(:, 2:end) ./ ENEL(:, 1:end-1);
ratio_AXA = AXA(:, 2:end) ./ AXA(:, 1:end-1);

% compute the final coupon
S_t = 1/4 * sum(1/2 * ratio_ENEL + 1/2 * ratio_AXA, 2);

% compute the discounted payoffs
disc_payoff = discounts_fixed(end) * max(S_t - P, 0) * principal_Amount;
% compute the mean
maturity_payment_B = mean(disc_payoff);
% initial party B payment
start_payment_B = X * principal_Amount;

%% Party A leg

% compute discount factors for the floating leg
discounts_floating = intExtDF(discounts, dates, floatingLegDates);

% compute the deltas
ACT_360 = 2;
deltas_floating = yearfrac(floatingLegDates(1:end-1), floatingLegDates(2:end), ACT_360);

% compute the s_spol payments
NPV_s_spol = s_spol * sum(deltas_floating .* discounts_floating(2:end)) * principal_Amount;

% compute the Libor payments
NPV_libor = (1 - discounts_fixed(end)) * principal_Amount;

% payment at maturity
maturity_payment_A = (1-P) * discounts_fixed(end) * principal_Amount;

party_A_leg = NPV_libor + NPV_s_spol + maturity_payment_A;

%% Compute the partecipation coefficient

alpha = (party_A_leg - start_payment_B) / maturity_payment_B

%% Point 2: Pricing Digital Option

% Price with Black Formula
Notional = 1e7;
% Underlying price
S_digital = cSelect.reference;
% Strike
k = S_digital;
% maturity
T = yearfrac(t0, t0 + calyears(1), ACT_365);
% value of the digital option at maturity if s > k
payoff = 0.05 * S_digital;
% compute the discount factor at 1 year
discount_1y = intExtDF(discounts, dates, datenum(t0 + calyears(1))); 

% Load volatility smile
strikes = cSelect.strikes;
surface = cSelect.surface;

sigma_digial = interp1(strikes, surface, k, 'spline');

plot(strikes, surface); hold on;
plot(k, sigma_digial,'x')

d_1 = (log(S_digital / k) + (S_digital + 0.5 * sigma_digial) * T) / (sigma_digial * sqrt(T));
d_2 = d_1 - sigma_digial * sqrt(T);

price_digital_black = payoff * discount_1y * normcdf(d_2);

price = price_digital_black * Notional


%% Point 3: Pricing


% parameters of the mean-variance mixture model
alpha = 0.5;
sigma = 20 / 100;
kappa = 1;
eta = 3;
t = 1;
% moneyness
x = -25:1:25 / 100;
%F_0 = 


%% Point 4: Volatility Surface Calibration

toc