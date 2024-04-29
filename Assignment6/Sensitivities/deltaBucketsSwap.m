function [sensitivities, sens_dates] = deltaBucketsSwap(datesSet, ratesSet, quoted_dates, swapRate, swapTtm)
% DELTABUCKETSSWAP Compute the delta-bucket sensitivities of a swap
%
% INPUTS
%   datesSet: dates of the market data
%   ratesSet: rates of the market data
%   quoted_dates: dates of the quoted rates to be shifted
%   swapRate: fixed rate of the swap
%   swapTtm: time to maturity of the swap

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(quoted_dates), 1);
sens_dates = quoted_dates;

% shock the mid-market rates by one basis point each and compute the
% change in NPV

% shift is 1 bp
shift = 0.0001; % 1 bp

% skip the first date (t0)
for i = 1:length(quoted_dates)
    % shift the corresponding rate by 1 bp
    shifted_ratesSet = shift_rate(ratesSet, datesSet, quoted_dates(i), shift);
    % rerun the bootstrap
    [shifted_dates, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    % compute the sensitivity
    sensitivities(i) = swapNPV(swapRate, swapTtm, shifted_discounts, shifted_dates);
end

end