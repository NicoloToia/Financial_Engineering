function [sens_dates, sensitivities] = deltaBucketsCap(cap_0, datesSet, ratesSet, quoted_dates, strike, ...
    ttm, spot_vols, spot_ttms, strikes)
% DELTABUCKETSCAP computes the delta-bucket sensitivities for a cap
%
% INPUTS
%   cap_0: initial price of the cap
%   datesSet: dates of the market data
%   ratesSet: rates of the market data
%   quoted_dates: dates of the quoted rates to be shifted
%   strike: strike of the cap
%   ttm: time to maturity of the cap
%   spot_vols: spot volatilities (calibrated from the market data)
%   spot_ttms: time to maturities of the spot_vols
%   strikes: strikes

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(quoted_dates), 1);
sens_dates = quoted_dates;

% compute the dates for the cap
t0 = datesSet.settlement;

exercise_dates = datetime(t0, 'ConvertFrom', 'datenum') + ...
    calmonths(3:3:ttm*12-3)';
payment_dates = exercise_dates + calmonths(3);
% move to business days if needed
exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
    busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
    busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
% conver to datenum
exercise_dates = datenum(exercise_dates);
payment_dates = datenum(payment_dates);

% shock the mid-market rates by one basis point each and compute the
% change in NPV

shift = 10^(-4); % 1 bp

% skip the first date (t0)
for i = 1:length(quoted_dates)
    % shift the rates of 1 bp in the bucket date
    shifted_ratesSet = shift_rate(ratesSet, datesSet, quoted_dates(i), shift);
    % rerun the bootstrap
    [shifted_dates, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    % compute the cap price for given maturity and strike
    cap_shift = CapSpot(strike, exercise_dates, payment_dates, spot_vols, spot_ttms, strikes, ...
        shifted_discounts, shifted_dates);
    % compute the delta-bucket sensitivity
    sensitivities(i) = (cap_shift - cap_0);
end

end