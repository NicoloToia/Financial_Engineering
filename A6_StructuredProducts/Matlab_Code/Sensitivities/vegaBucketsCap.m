function sensitivities = vegaBucketsCap(cap_price_0, strike, ttm, mkt_vols, ttms, strikes, discounts, dates)
% VEGABUCKETSCAP computes the vega bucket sensitivities for the certificate
%
% INPUTS
%   cap_price_0: initial price of the cap
%   strike: strike of the cap
%   ttm: time to maturity of the cap
%   mkt_vols: market volatilities
%   ttms: time to maturities of the mkt_vols
%   strikes: strikes
%   discounts: discounts
%   dates: dates of the market data

% initialize the sensitivities
sensitivities = zeros(length(mkt_vols), 1);
% shift is 1 bp
shift = 10^(-4);

% compute the cap dates
exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:ttm*12-3)';
payment_dates = exercise_dates + calmonths(3);
% move to business days if needed
exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
    busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
    busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
% convert to datenum
exercise_dates = datenum(exercise_dates);
payment_dates = datenum(payment_dates);

% for each maturity, compute the vega bucket sensitivity
for i = 1:length(mkt_vols)
    % shift the row by 1 bp
    [shifted_ttms, shifted_vols] = shiftVolsRow(mkt_vols, i, shift, ttms, strikes, discounts, dates);
    % recompute the cap price
    Cap_shift = CapSpot(strike, exercise_dates, payment_dates, shifted_vols, shifted_ttms, strikes, discounts, dates);
    % compute the vega bucket sensitivity
    sensitivities(i) = (Cap_shift - cap_price_0);
end

end