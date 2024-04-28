function X = computeUpfront(spot_vols, ttms, strikes, start_date, spol_A, fixed_rate_B,...
    spol_B, cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates)
% computeUpfront: Compute the upfront of the certificate (see annex for more details)
%
% Annex
% -----
% Maturity: 15 years
%
% - Party A payments:
%   - Euribor 3m + 2.00% quarterly daycount ACT/360
% - Party B payments:
%   - X% of principal at start date
%   - 3.00% on the first quarter
%   - On the following quarters
%       - min(Euribor 3m + 1.10%, 4.30%) up to and including 5 years
%       - min(Euribor 3m + 1.10%, 4.60%) from 5 years up to and including 10 years
%       - min(Euribor 3m + 1.10%, 5.10%) from 10 years up to and including 15 years
%
% INPUT
%   spot_vols : Spot volatilities
%   ttms : Time to maturities of the caplets
%   strikes : Strikes of the caplets
%   start_date : Start date of the certificate
%   spol_A : flat rate that party A pays on top of the Euribor 3m
%   fixed_rate_B : fixed rate that party B pays on the first quarter of each year
%   spol_B : flat rate that party B pays on top of the Euribor 3m
%   cap_rate_5y : cap rate for the first 5 years
%   cap_rate_10y : cap rate for the first 10 years
%   cap_rate_15y : cap rate for the first 15 years
%   discounts : Discount factors curve
%   dates : Dates of the discount factors

% set the day count convention
ACT_360 = 2;

% compute the payment dates and exercise dates

% exercise dates
exercise_dates = datetime(start_date, 'ConvertFrom', 'datenum') + calmonths(0:3:12*15-3)';
payments_dates = exercise_dates + calmonths(3)';

% move to business days
exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
    busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
payments_dates(~isbusday(payments_dates, eurCalendar())) = ...
    busdate(payments_dates(~isbusday(payments_dates, eurCalendar())), 'modifiedfollow', eurCalendar());

% convert to datenum
exercise_dates = datenum(exercise_dates);
payments_dates = datenum(payments_dates);

% compute the year fractions between the exercise dates and payment dates
deltas_libor = yearfrac(exercise_dates, payments_dates, ACT_360);

% compute the discount factors at the payment dates and exercise dates
DF_payment = intExtDF(discounts, dates, payments_dates);

% party A payments
Libor_payment_A = 1 - DF_payment(end);
spol_payment_A = spol_A / 100 * deltas_libor' * DF_payment;

NPV_A = Libor_payment_A + spol_payment_A;

% party B payments
% the following relation holds: min(x+q, k) = x + q - max(0, x+q-k) = x + q - caplet(x, k-q)
% our payments factorize
% - First term:
%   this is simply the libor for each date from the t0+3m onwards
%   hence the libor payments simplify to: 1 - B(0,15y) - B(0,3m) * delta(0,3m) * L(0,3m)
% - Second term:
%   this is just the spol of B times the bpv from 6m to 15y quarter by quarter
% - Third term:
%   this is the caplet payments
%   - From 0 to 5y: The caplets sum and form the cap with strike K = 4.3% - 1.1% from 0 to 5y
%   - From 5 to 10y: The caplets sum and form the cap with strike K = 4.6% - 1.1% from 5 to 10y
%       this can further be written as the difference of the 0 to 10y cap minus the 0 to 5y cap
%       with the same strike
%   - From 10 to 15y: The caplets sum and form the cap with strike K = 5.1% - 1
%       this can further be written as the difference of the 0 to 15y cap minus the 0 to 10y cap

% first quarter fixed rate payment
first_quarter_B = fixed_rate_B / 100 * deltas_libor(1) * DF_payment(1);

% Libor payments
Libor_payment_B = DF_payment(1) - DF_payment(end);

% fixed rate payments
fixed_rate_payment_B = spol_B / 100 * deltas_libor(2:end)' * DF_payment(2:end);

% caplet payments

% from 0 to 5y
% compute the strike
strike_5y = cap_rate_5y - spol_B;
% find the relevant dates (from 2nd quarter to 5y)
exercise_dates_5y = exercise_dates(2:4*5);
payments_dates_5y = payments_dates(2:4*5);

% compute the cap from 0 to 5y with given strike
cap_5y = CapSpot(strike_5y, exercise_dates_5y, payments_dates_5y, spot_vols, ttms, strikes, ...
    discounts, dates);

% from 5 to 10y
% compute the strike
strike_10y = cap_rate_10y - spol_B;
% find the relevant dates (from 5y to 10y)
exercise_dates_10y = exercise_dates(4*5+1:4*10);
payments_dates_10y = payments_dates(4*5+1:4*10);

% compute the cap from 5 to 10y with given strike
cap_10y = CapSpot(strike_10y, exercise_dates_10y, payments_dates_10y, spot_vols, ttms, strikes, ...
    discounts, dates);

% from 10 to 15y
% compute the strike
strike_15y = cap_rate_15y - spol_B;
% find the relevant dates (from 10y to 15y)
exercise_dates_15y = exercise_dates(4*10+1:end);
payments_dates_15y = payments_dates(4*10+1:end);

% compute the cap from 10 to 15y with given strike
cap_15y = CapSpot(strike_15y, exercise_dates_15y, payments_dates_15y, spot_vols, ttms, strikes, ...
    discounts, dates);

% compute the upfront
X = NPV_A - (first_quarter_B + Libor_payment_B + fixed_rate_payment_B - cap_5y - cap_10y - cap_15y);

end