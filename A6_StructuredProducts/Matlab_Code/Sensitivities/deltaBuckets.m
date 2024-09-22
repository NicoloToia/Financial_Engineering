function [sens_dates, sensitivities] = deltaBuckets(datesSet, ratesSet, quoted_dates, spot_vols, ttms, strikes, X_0, ...
    start_date, spol_A, fixed_rate_B, spol_B, cap_5y, cap_10y, cap_15y)
% DELTABUCKETS computes the delta-bucket sensitivities
%
% INPUTS
%   datesSet: dates of the market data
%   ratesSet: rates of the market data
%   quoted_dates: dates of the quoted rates to be shifted
%   spot_vols: spot volatilities
%   ttms: time to maturities
%   strikes: strikes
%   X_0: upfront payment in the initial market data
%   start_date: start date of the market data
%   spol_A: spread over libor for the first quarter
%   fixed_rate_B: fixed rate for the second quarter
%   spol_B: spread over libor for the third quarter
%   cap_5y: strike of the cap from 0 to 5y
%   cap_10y: strike of the cap from 5 to 10y
%   cap_15y: strike of the cap from 10 to 15y

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(quoted_dates), 1);
sens_dates = quoted_dates;

% shock the mid-market rates by one basis point each and compute the
% change in NPV

% shift is 1 bp
shift = 10^(-4); % 1 bp

% skip the first date (t0)
for i = 1:length(quoted_dates)
    % shift the rates of 1 bp in the bucket date
    shifted_ratesSet = shift_rate(ratesSet, datesSet, quoted_dates(i), shift);
    % rerun the bootstrap
    [shifted_dates, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    X_shift = computeUpfront(spot_vols, ttms, strikes, start_date, spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, shifted_discounts, shifted_dates);
    % compute the delta-bucket sensitivity
    sensitivities(i) = (X_shift - X_0);
end

end