function [X_hat, IC] = price3y(S_0, d, K, ttm, alpha, sigma, kappa, eta, s_A, N, discounts, dates, ...
    principal, coupons, alpha_IC)
% price3y: function to price a 3y autocallable certificate
% INPUTS
%   S0: initial price of the underlying
%   d: dividend yield rate
%   K: strike price
%   ttm: time to maturity of the certificate
%   alpha = parameter to identify the desired model, in this 0.5(NIG)
%   sigma: volatility of the underlying
%   kappa: vol of vol
%   eta: skewness
%   s_A: spread for party A
%   N = number of simulations for MC
%   discounts: discount factors
%   dates: dates for the discount factors
%   principal : principal amount
%   coupons : vector of possible coupons
%   alpha_IC : confidence level for the IC

% Define yearfractions conventions
ACT_360 = 2; 
ACT_365 = 3; % Actual/365 day count convention
EU_30_360 = 6; % 30/360 day count convention

% compute the coupon dates (one each year)
couponDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(1:ttm)';
couponDates(~isbusday(couponDates, eurCalendar())) = ...
    busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);

% compute the yearfracs and discount factors at coupon dates
deltas_coupon = yearfrac([dates(1); couponDates(1:end-1)], couponDates, EU_30_360);
coupon_DF = intExtDF(discounts, dates, couponDates);

% find the date up to which to simulate the underlying and the discount factors
fixing_dates = datetime(couponDates, 'ConvertFrom', 'datenum') - caldays(2);
fixing_dates(~isbusday(fixing_dates, eurCalendar())) = ...
    busdate(fixing_dates(~isbusday(fixing_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
fixing_dates = datenum(fixing_dates);
fixing_DF = intExtDF(discounts, dates, fixing_dates);

% compute the dates, yearfracs and discount factors for the quarterly payments
quarter_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12*ttm)';
quarter_dates(~isbusday(quarter_dates, eurCalendar())) = ...
    busdate(quarter_dates(~isbusday(quarter_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
quarter_dates = datenum(quarter_dates);
quarter_DF = intExtDF(discounts, dates, quarter_dates);
quarter_deltas = yearfrac([dates(1); quarter_dates(1:end-1)], quarter_dates, ACT_360);

% simulation
t = yearfrac(dates(1), fixing_dates, ACT_365); % find the time step for the simulation
S = zeros(N, length(t)); % initialize the values of the underlying
F = S_0 / fixing_DF(end) * exp(-d * t(end)); % initial price for the forward 

% simulate the forward and compute the underlying price
for i = 1:length(t)
    if i == 1
        dt = t(i);
    else
        dt = t(i) - t(i-1);
    end
    % simulate the forward
    F = MC_NIG(F, alpha, sigma, kappa, eta, dt, N);
    % compute the value of the underlying
    S(:, i) = F * fixing_DF(end) / fixing_DF(i) * exp(d * (t(end) - t(i)));
end

% compute the prices for the coupons
first_coupon = coupons(1) * deltas_coupon(1) * coupon_DF(1) * (S(:, 1) < K);
second_coupon = coupons(2) *  deltas_coupon(2) * coupon_DF(2) * (S(:, 2) < K);
third_coupon = coupons(3) *  deltas_coupon(3) * coupon_DF(3) * (S(:, 2) > K) .* (S(:, 1) > K );

% NPV_B 
NPV_B = principal * (first_coupon + second_coupon + third_coupon);

% compute the party A leg

% compute the BPV for 1 and 2 years
BPV_2y = sum(quarter_deltas(1:8) .* quarter_DF(1:8));
BPV_3Y = sum(quarter_deltas(9:12) .* quarter_DF(9:12)); 

% the payments for the the first two years get paid with certainty.
% the last only if the third coupon is paid
NPV_A = principal * (1 - quarter_DF(8) + s_A * BPV_2y) + ...
    principal * (quarter_DF(8) - quarter_DF(12) + s_A * BPV_3Y) * (S(:, 2) > K) .* (S(:, 1) > K );

% compute the upfront for each simulation
X = NPV_A - NPV_B;

% take the mean and variance to compute the IC
X_hat = mean(X);
sigma = std(X);

% confidence interval at 95%
z = norminv(1- alpha_IC / 2);
IC = [X_hat - z * sigma / sqrt(N), X_hat + z * sigma / sqrt(N)];

end
