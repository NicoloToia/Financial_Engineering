function X = computeUpfrontSkew(S0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, ...
    discounts, dates, strikes, surface, flag)
% computeUpfrontSkew: computes the upfront payment for a certificate with a skew
%
% INPUTS
%   S0: initial price of the underlying
%   d: dividend yield rate
%   K: strike price
%   ttm: time to maturity of the certificate
%   coupon_1y: coupon rate for the 1st year payment
%   coupon_2y: coupon rate for the 2nd year payment
%   s_A: spread for party A
%   sigma: volatility of the underlying
%   kappa: vol of vol
%   eta: skewness
%   discounts: discount factors
%   dates: dates for the discount factors

% Define yearfractions conventions
ACT_360 = 2; % Actual/360 day count convention
ACT_365 = 3; % Actual/365 day count convention
EU_30_360 = 6; % 30/360 day count convention

% compute the coupon dates
couponDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(1:ttm)';
couponDates(~isbusday(couponDates, eurCalendar())) = ...
    busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);

% compute the deltas
delta_1y = yearfrac(dates(1), couponDates(1), EU_30_360);
delta_2y = yearfrac(couponDates(1), couponDates(2), EU_30_360);

% compute the discount factors at the payment dates
coupon_DF = intExtDF(discounts, dates, couponDates);

% find the date up to which to simulate the underlying
fixing_date = datetime(couponDates(1), 'ConvertFrom', 'datenum') - caldays(2);
fixing_DF = intExtDF(discounts, dates, datenum(fixing_date));
t = yearfrac(dates(1), fixing_date, ACT_365);

% initial data
F_0 = S0 / fixing_DF * exp(-d * t);

% compute the digital call price
price_digital = DigitalPrice(1, t, F_0, d, fixing_DF, strike, strikes, surface, flag);

% compute the probability that the underlying is above the strike
prob_up = price_digital / fixing_DF;

% compute the party B leg
NPV_B = principal * delta_1y * coupon_1y * coupon_DF(1) * (1 - prob_up) + ...
    principal * delta_2y * coupon_2y * coupon_DF(2) * prob_up;

% compute the party A leg

% compute the quarterly payments
quarter_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12*ttm)';
quarter_dates(~isbusday(quarter_dates, eurCalendar())) = ...
    busdate(quarter_dates(~isbusday(quarter_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
quarter_dates = datenum(quarter_dates);

% compute the discount factors at the payment dates
quarter_DF = intExtDF(discounts, dates, quarter_dates);
% compute the deltas for the BPV
deltas = yearfrac([dates(1); quarter_dates(1:end-1)], quarter_dates, ACT_360);

% compute the BPV for 1 and 2 years
BPV_1y = sum(deltas(1:4) .* quarter_DF(1:4));
BPV_2y = sum(deltas(5:8) .* quarter_DF(5:8));

% compute the NPV for party A
NPV_A = principal * (1 - quarter_DF(4) + s_A * BPV_1y) + ...
    principal * (quarter_DF(4) - quarter_DF(8) + s_A * BPV_2y) * prob_up;

% compute the upfront
X = NPV_A - NPV_B;

end