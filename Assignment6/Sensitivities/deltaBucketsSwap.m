function [sensitivities, sens_dates] = deltaBucketsSwap(datesSet, ratesSet, dates, swapRate, swapDates)
% DELTABUCKETSSWAP Compute the delta-bucket sensitivities of a swap
%
% INPUTS
%   datesSet: dates of the market data
%   ratesSet: rates of the market data
%   dates: dates to shock
%   swapRate: fixed rate of the swap
%   swapDates: dates of the swap (including the settlement date)

% initialize the sensitivities (skip the first date) and the dates
sensitivities = zeros(length(dates), 1);
sens_dates = dates;

% shock the mid-market rates by one basis point each and compute the
% change in NPV

shift = 0.0001; % 1 bp

% skip the first date (t0)
for i = 2:length(dates)
    % shift the rates of 1 bp in the bucket date
    try
        shifted_ratesSet = shift_rate(ratesSet, datesSet, dates(i), shift);
    catch
        % remove the date from the sensitivities and carry on
        sensitivities(i) = NaN;
        continue;
    end
    % rerun the bootstrap
    [~, shifted_discounts] = bootstrap(datesSet, shifted_ratesSet);
    sensitivities(i) = swapNPV(swapRate, swapDates, shifted_discounts, dates) / shift;
end

% remove the NaNs and the corresponding dates
sens_dates = sens_dates(~isnan(sensitivities));
sensitivities = sensitivities(~isnan(sensitivities));

% remove the zero sensitivities
% sens_dates = sens_dates(sensitivities ~= 0);
% sensitivities = sensitivities(sensitivities ~= 0);

end