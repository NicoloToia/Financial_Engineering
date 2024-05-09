function X = price3y(S_0, d, K, ttm, alpha, sigma, kappa, eta, s_A, N, discounts, dates, principal, coupons)

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
%   principal = principal amount
%   coupons = vector of possible coupons


% Define yearfractions conventions
ACT_360 = 2; 
ACT_365 = 3; % Actual/365 day count convention
EU_30_360 = 6; % 30/360 day count convention

% compute the coupon dates
couponDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(1:ttm)';
couponDates(~isbusday(couponDates, eurCalendar())) = ...
    busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);

% compute the deltas
deltas(1) = yearfrac(dates(1), couponDates(1), EU_30_360);
deltas(2) = yearfrac(couponDates(1), couponDates(2), EU_30_360);
deltas(3) = yearfrac(couponDates(2), couponDates(3), EU_30_360);

% compute the discount factors at the payment dates
coupon_DF = intExtDF(discounts, dates, couponDates);

% find the date up to which to simulate the underlying
fixing_date(1) = datetime(couponDates(1), 'ConvertFrom', 'datenum') - caldays(2);
fixing_date(2) = datetime(couponDates(2), 'ConvertFrom', 'datenum') - caldays(2);

%%%%%%% business day check?

%find the discount factors via linear interpolation
fixing_DF = intExtDF(discounts, dates, datenum(fixing_date));

%find the yearfrac between the fixing dates
t(1) = yearfrac(dates(1), fixing_date(1), ACT_365);
t(2) = yearfrac(fixing_date(1), fixing_date(2), ACT_365);

% initial data
F_0 = S_0 / fixing_DF(2) * exp(-d * (t(1) + t(2)) );

% compute the forward price via Monte Carlo NIG simulation
FT_1 = MC_NIG(F_0, alpha, sigma, kappa, eta, t(1), N);

% Compute the undelying price
ST_1 = FT_1* fixing_DF(1) * exp(d*t(1));

% Divide ST in two vectors: one with ST>K and the other one with ST<K
FT_1up = FT_1(ST_1>K);
%I don't take the first coupon
FT_1down = FT_1(ST_1<K); %I take the first coupon

%Compute the forward price via MonteCarlo NIG starting from FT_1
FT_2up = MC_NIG(FT_1up, alpha, sigma, kappa, eta, t(2), length(FT_1up)); 
FT_2down = MC_NIG(FT_1down, alpha, sigma, kappa, eta, t(2), length(FT_1down)); 

% Compute the undelying price
ST_2up = FT_2up;
ST_2down = FT_2down;

%now we find 4 different scenarios: 
%for each scenario, find the probability and the corresponding payoff

%case A: ST_1>K & ST_2>K -> we take the last coupon
ST_A = ST_2up(ST_2up>K); 
probA = length(ST_A)/N; 
partA = principal * deltas(3) * coupons(3) * coupon_DF(3) * probA;

%case B: ST_1>K & ST_2<K -> we take the second coupon
ST_B = ST_2up(ST_2up<K); 
probB = length(ST_B)/N; 
partB = principal * deltas(2) * coupons(2) * coupon_DF(2) * probB; 

%case C: ST_1<K & ST_2>K -> we take the first 
ST_C = ST_2down(ST_2down>K); 
probC = length(ST_C)/N; 
partC = principal * deltas(1) * coupons(1) * coupon_DF(1) * probC; 

%case D: ST_1<K & ST<K -> we take the first and the second coupon 
ST_D = ST_2down(ST_2down<K); 
probD = length(ST_D)/N; 
partD = principal * (deltas(1) * coupons(1) * coupon_DF(1) + ...
    deltas(2) * coupons (2) * coupon_DF(2)) * probD; 

%compute the prices
%NPV_B 
NPV_B = partA + partB + partC + partD;  

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
BPV_3Y = sum(deltas(9:12) .* quarter_DF(9:12)); 

% compute the NPV for party A
NPV_A = principal * (1 - quarter_DF(4) + s_A * BPV_1y) + ...
    principal * (quarter_DF(4) - quarter_DF(8) + s_A * BPV_2y) + ...
    principal * (quarter_DF(8) - quarter_DF(12) + s_A * BPV_3Y) * (probA); 

% compute the upfront
X = NPV_A - NPV_B;

end





