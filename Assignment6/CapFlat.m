function Price = CapFlat(Strike, Vol, startDate, paymentDates, skipFirst, discounts, dates)
% CapFlat: Compute the price of a cap using the flat volatility assumption
%
% INPUT
%   Strike : Strike price of the cap
%   Vol : Volatility of the cap
%   paymentDates : Payment dates of the cap (do not include the settlement date)
%   skipFirst : Skip the first payment date
%   discounts : Discount factors
%   dates : Dates of the discount factors

% compute the deltas
ACT_360 = 2;
deltas = yearfrac([startDate;paymentDates(1:end-1)], paymentDates, ACT_360);

% compute the discount factors
DF = intExtDF(discounts, dates, [startDate;paymentDates]); % this is one more than the deltas

% compute the forward rates
fwdDiscounts = DF(2:end) ./ DF(1:end-1);
% compute the forward Libor
% the ith represents the forward Libor from i to i+1
fwdLibor = 1 ./ deltas .* (1 ./ fwdDiscounts - 1);

% transform the skipFirst from boolean to integer
if skipFirst
    skipFirst = 1;
else
    skipFirst = 0;
end

% compute the caplets prices using the Bachelier formula
Caplets = zeros(length(paymentDates), 1);

% price the first Caplet
for i = 1+skipFirst:length(paymentDates)

    Caplets(i) = CapletBachelier(fwdLibor(i), (1+Strike/100)*fwdLibor(i), Vol, deltas(i), paymentDates(i), dates(1), DF(i+1));
    
end

Price = sum(Caplets);

end