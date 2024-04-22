function CapPrice = Cap3mFlat(Strike, Vol, paymentDates, startDate, skipFirst, discounts, dates)

% compute the deltas
ACT_360 = 2;
deltas = yearfrac(paymentDates(1:end-1), paymentDates(2:end), ACT_360);

% compute the discount factors
DF = intExtDF(discounts, dates, [startDate; paymentDates]);
% if any NaN is present in the discount factors, substitute with 1
DF(isnan(DF)) = 1;

% compute the forward rates
fwdDiscounts = DF(2:end) ./ DF(1:end-1);
% compute the forward Libor
fwdLibor = 1 ./ deltas .* (1 ./ fwdDiscounts - 1);

% skip the first PaymentDate if necessary
if skipFirst
    paymentDates = paymentDates(2:end);
    deltas = deltas(2:end);
    fwdLibor = fwdLibor(2:end);
    DF = DF(2:end);
end

% compute the caplets prices using the Bachelier formula
Caplets = zeros(length(paymentDates), 1);

% price the first Caplet
Caplets(1) = CapletBachelier(fwdLibor(1), Strike, Vol, deltas(1), startDate, dates(1), DF(1));
for i = 2:length(paymentDates)
    Caplets(i) = CapletBachelier(fwdLibor(i), Strike, Vol, deltas(i), paymentDates(i-1), dates(1), discounts(i));
end

CapPrice = sum(Caplets);

end