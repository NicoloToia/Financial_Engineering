function CapPrices = capPrices(ttms, strikes, Vols, discounts, dates)
%

% save the caps
CapPrices = zeros(length(ttms), length(strikes));

for j = 1:length(ttms)
    for i = 1:length(strikes)

        first_year_payment_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:(12*j))';
        first_year_payment_dates = datenum(first_year_payment_dates);
        CapPrices(j,i) = CapFlat(strikes(i), Vols(j, i), dates(1), first_year_payment_dates, true, discounts, dates);

    end
    
end

end