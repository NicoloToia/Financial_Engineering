function [discounts] = discountSwaps(datesSet, ratesSet, discounts, dates , SwapDate, swapStart)
% Bootstrap the discount factors for swaps
%
% INPUT
% datesSet   : dates from market data stored as follow
%              -> datesSet.settlement : Settlement Date
%              -> datesSet.depos      : Vector of end calculation dates for Depos
%              -> datesSet.futures    : Matrix with start/end dates for Futures
%              -> datesSet.swaps      : Vector of end calculation date for Swaps
% ratesSet   : market data for bid&ask rates stored as follow
%              -> ratesSet.depos      : Matix bid&ask Depos
%              -> ratesSet.futures    : Matix bid&ask Futures
%              -> ratesSet.swaps      : Matix bid&ask Swaps
% discounts  : Inizialized vector for discounts of the length needed, 
%             set as first value one
% dates      : dates of discount factors curve bootstrap
% SwapStart : the function compute the discounts from the n swap,
%             where n = SwapStart
%
% OUTPUT
% discounts: Discount Factors for swaps


% Yearfrac Convention 
EU_30_360 = 6;

% Retrieve rates and dates
swapsDates = datesSet.swaps;                            % first 7 swaps dates
swapsRates = ratesSet.swaps;                            % first 7 swaps rates

% Fix settlement date
t0 = datesSet.settlement;

% Interpolate the first discount factor (swap 1y) and save it
delta = yearfrac(t0, swapsDates(1), EU_30_360);
% use only previous values
discFact = intExtDF(discounts(1:swapStart-1), dates(1:swapStart-1), swapsDates(1));

% BPV represents BPV(0, T_i-1)
% Initialize with the first Discount factor
BPV = delta * discFact;

% compute the discount factors for the swaps
for i=SwapDate:length(swapsDates)
    % compute delta for t_i-1, t_i
    delta = yearfrac(swapsDates(i-1), swapsDates(i), EU_30_360);
    % compute new discount factor for t_i
    discount = (1 - swapsRates(i) * BPV) / (1 + swapsRates(i) * delta);
    % update discounts
    discounts(swapStart + i - 2) = discount;
    % update the BPV
    BPV = BPV + delta * discount;
end

end