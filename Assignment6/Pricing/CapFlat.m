function Price = CapFlat(Strike, Vol, exercise_dates, payments_dates, skipFirst, discounts, dates)
% CapFlat: Compute the price of a cap using the flat volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%)
%   Vol : Volatility of the cap
%   exercise_dates : Dates in which the caplets are exercised
%   payment_dates : Dates in which the Euribor is paid
%   skipFirst : Skip the first caplet
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

% transform the skipFirst from boolean to integer
if skipFirst
    skipFirst = 1;
else
    skipFirst = 0;
end

% compute the caplets prices using the Bachelier formula
Caplets = zeros(length(payments_dates), 1);

% price the first Caplet
for i = 1+skipFirst:length(payments_dates)

    % compute the price of the caplet
    Caplets(i) = CapletBachelier(fwd_libor(i), Strike, Vol, payments_dates(i), exercise_dates(i), ...
        dates(1), DF_payment(i));
    
end

% sum the caplets to get the cap price
Price = sum(Caplets);

end