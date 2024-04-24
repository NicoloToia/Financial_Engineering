function Price = CapSpot(strike, exercise_dates, payment_dates, spot_vols, ttms, strikes, ...
    discounts, dates)
% CapSpot: Compute the price of a cap using the spot volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%) 
%   exercise_dates : Dates in which the caplets are exercised
%   payment_dates : Dates in which the Euribor is paid
%   spot_vols : Spot volatilities (calibrated to the market caplets)
%   ttms : Time to maturities of the caplets (in years)
%   strikes : Strikes of the caplets
%   discounts : Discount factors curve
%   dates : Dates of the discount factors

% find the indices of the two columns of the spot volatilities that are closest to the
% target strike
idx1 = find(strikes <= strike, 1, 'last');
idx2 = find(strikes >= strike, 1, 'first');

strike_1 = strikes(idx1)
strike_2 = strikes(idx2)


end