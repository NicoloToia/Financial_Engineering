function [Price, sigmas] = CapSpotBootStrap(Strike, sigma_alpha, T_alpha, sigma_beta, exercise_dates, payments_dates, discounts, dates)
% CapFlat: Compute the price of a cap using the flat volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%)
%   sigma_alpha : Volatility of the previous caplet
%   sigma_beta : Volatility of final caplet
%   exercise_dates : Dates in which the caplets are exercised
%   payment_dates : Dates in which the Euribor is paid
%   discounts : Discount factors curve
%   dates : Dates of the discount factors

% compute the discount factor at the payment dates
DF_payment = intExtDF(discounts, dates, payments_dates);

% compute the Libor
fwd_discounts = DF_payment ./ intExtDF(discounts, dates, exercise_dates);

% compute the forward Libor
% the ith represents the forward Libor from i to i+1
ACT_360 = 2;
deltas = yearfrac(exercise_dates, payments_dates, ACT_360);
fwd_libor = 1 ./ deltas .* (1 ./ fwd_discounts - 1);

ACT_365 = 3;
sigmas = sigma_alpha + (sigma_beta - sigma_alpha) / yearfrac(T_alpha, exercise_dates(end), ACT_365) * ...
    yearfrac(T_alpha, exercise_dates, ACT_365);

% compute the caplets prices using the Bachelier formula
Caplets = zeros(length(payments_dates), 1);

% price the first Caplet
for i = 1:length(payments_dates)

    Caplets(i) = CapletBachelier(fwd_libor(i), Strike, sigmas(i), payments_dates(i), exercise_dates(i), ...
        dates(1), DF_payment(i));
    
end

Price = sum(Caplets);

end