function CapPrices = capPrices(ttms, strikes, Vols, discounts, dates)
%

% save the caps
CapPrices = zeros(length(ttms), length(strikes));

% for each strike price the first cap
first_year_payment_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12)';
first_year_payment_dates = datenum(first_year_payment_dates);
for i = 1:length(strikes)
    CapPrices(1,i) = Cap3mFlat(strikes(i), Vols(1, i), first_year_payment_dates, dates(1), true, discounts, dates);
end

CapPrices

% price the caps
for i = 2:length(ttms)
    for j = 1:length(strikes)
    end
end
    
end