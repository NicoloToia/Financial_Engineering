function [sens_dates, sensitivities] = deltaBucketsCap(cap_0, datesSet, ratesSet, quoted_dates, strike, ...
    ttm,  caplet_ttms, caplet_yf, exercise_dates, payment_dates, spot_vols, spot_ttms, strikes)
% DELTABUCKETSCAP computes the delta-bucket sensitivities for a cap
%
% INPUTS
%   cap_0: initial price of the cap
%   datesSet: dates of the market data
%   ratesSet: rates of the market data
%   quoted_dates: dates of the quoted rates to be shifted
%   strike: strike of the cap
%   ttm: time to maturity of the cap
%   caplet_ttms: caplet maturities
%   caplet_yf: caplet year fractions
%   exercise_dates: exercise dates of the caplets
%   payment_dates: payment dates of the caplets
%   spot_vols: spot volatilities (calibrated from the market data)
%   spot_ttms: time to maturities of the spot_vols
%   strikes: strikes

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(quoted_dates), 1);
sens_dates = quoted_dates;

% find the relevant ttms and yf for the caplets
cap_caplet_yf = caplet_yf(2:15*4);
cap_caplet_ttms = caplet_ttms(2:15*4);

% shock the mid-market rates by one basis point each and compute the
% change in NPV

shift = 10^(-4); % 1 bp

% skip the first date (t0)
for i = 1:length(quoted_dates)
    % shift the rates of 1 bp in the bucket date
    shifted_ratesSet = shift_rate(ratesSet, datesSet, quoted_dates(i), shift);
    % rerun the bootstrap
    [shifted_dates, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    % compute the caplet discount factors and the Libor rates
    caplet_DF = intExtDF(shifted_discounts, shifted_dates, payment_dates);
    fwd_DF = caplet_DF ./ intExtDF(shifted_discounts, shifted_dates, exercise_dates);
    Libor = (1 ./ fwd_DF - 1) ./ caplet_yf;
    % use only the relevant caplet data
    cap_Libor = Libor(2:15*4);
    cap_caplet_DF = caplet_DF(2:15*4);
    % compute the cap price for given maturity and strike
    cap_shift = CapSpot(strike, cap_caplet_ttms, cap_caplet_yf, cap_caplet_DF, cap_Libor, spot_vols, spot_ttms, strikes);
    % compute the delta-bucket sensitivity
    sensitivities(i) = (cap_shift - cap_0);
end

end