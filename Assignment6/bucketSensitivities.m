function sens = bucketSensitivities(datesSet, ratesSet, dates, discounts)

% find the dates for the computation of the bucket sensitivities
bucketDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears([2, 5, 10, 15, 20, 30, 40, 50])';

bp = 0.0001;

% shift the the zero rates curve by 1bp on the bucket dates


% Compute Delta-bucket sensitivities

for i = 1 : length(bucketDates)
    % shift the rates
    rates_Set_shift = shift_rates(ratesSet, bucketDates(i),bp);
    % Bootstrap the curve with the shifted rates
    [dates_shift, discounts_shift] = bootstrap(datesSet, rates_Set_shift, datesSet);
    % Compute the bucket sensitivities
    sens(i) = abs();
end

end
