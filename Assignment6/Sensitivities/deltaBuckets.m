function [sens_dates, sensitivities] = deltaBuckets(datesSet, ratesSet, quoted_dates, spot_vols, ttms, strikes, X_0, ...
    spol_A, fixed_rate_B, spol_B, cap_5y, cap_10y, cap_15y, caplet_ttms, caplet_yf, exercise_dates, payment_dates)
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
%   caplet_ttms: caplet maturities
%   caplet_yf: caplet year fractions
%   caplet_DF: caplet discount factors

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(quoted_dates), 1);
sens_dates = quoted_dates;

% shock the mid-market rates by one basis point each and compute the
% change in NPV

% shift is 1 bp
shift = 10^(-4); % 1 bp

% save the relevant caplet year fractions and ttms
cf_caplet_yf = caplet_yf(1:15*4);
cf_caplet_ttms = caplet_ttms(1:15*4);

% skip the first date (t0)
for i = 1:length(quoted_dates)
    % shift the rates of 1 bp in the bucket date
    shifted_ratesSet = shift_rate(ratesSet, datesSet, quoted_dates(i), shift);
    % rerun the bootstrap
    [shifted_dates, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    % recompute the discount factors and the libor rates
    shifted_caplet_DF = intExtDF(shifted_discounts, shifted_dates, payment_dates);
    shifted_fwd_DF = shifted_caplet_DF ./ intExtDF(shifted_discounts, shifted_dates, exercise_dates);
    shifted_Libor = (1 ./ shifted_fwd_DF - 1) ./ caplet_yf
    % use only the relevant caplet data
    shifted_Libor = shifted_Libor(1:15*4);
    shifted_caplet_DF = shifted_caplet_DF(1:15*4);
    % compute the upfront payment
    X_shift = computeUpfront(spot_vols, ttms, strikes, spol_A, fixed_rate_B, spol_B, cap_5y, cap_10y, cap_15y, ...
        cf_caplet_ttms, cf_caplet_yf, shifted_caplet_DF, shifted_Libor);
    % compute the delta-bucket sensitivity
    sensitivities(i) = (X_shift - X_0);
end

end