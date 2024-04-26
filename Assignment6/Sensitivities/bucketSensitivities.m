function [delta_bucket, vega_bucket] = bucketSensitivities(datesSet, ratesSet, dates, discounts, spot_vols, spol_A, fixed_rate_B, ...
    spol_B, cap_5y, cap_10y, cap_15y, ttms, strikes)

%bucket years for the computation of the bucket sensitivities
years_bucket = [2, 5, 10, 15, 20, 30, 40, 50];
% find the dates for the computation of the bucket sensitivities
bucketDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(years_bucket)';
bp = 0.0001;
delta_bucket = zeros (1, length(bucketDates));
vega_bucket = zeros (1, length(bucketDates));

for i = 1 : length(bucketDates)
    % Compute Delta-bucket sensitivities

    % shift the rates of 1 bp in the bucket date
    rates_Set_shift = shift_rates(ratesSet, datesSet, bucketDates(i),bp);
    % Bootstrap the curve with the shifted rates
    [dates_shift, discounts_shift] = bootstrap(datesSet, rates_Set_shift);
    %we have to compute the prices of the certificate with and without the shifted curve
    price_shift = computeUpfront(spot_vols, ttms, strikes, dates_shift(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts_shift, dates_shift);
    price_0 = computeUpfront(spot_vols, ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    % Compute the bucket sensitivities with the shifted curve
    delta_bucket(i) = abs(price_shift - price_0);

    % Compute the vega bucket 
    
    %shift the spot volatility of 1% in the bucket date
    vols_shift = shift_vols(spot_vols, datesSet, bucketDates(i), 0.01);
    %we have to compute the prices of the certificate with and without the shifted volatility
    price_shift = computeUpfront(vols_shift, ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    price_0 = computeUpfront(spot_vols, ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    % Compute the bucket sensitivities with the shifted volatility
    vega_bucket(i) = (price_shift - price_0)/price_0;
end
end
