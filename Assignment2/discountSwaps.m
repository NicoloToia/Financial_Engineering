function [discounts dates] = discountSwaps(datesSet, ratesSet, discounts, dates, SwapStart)

% date notation 
EU_30_360 = 6;

% retrieve data
swapsDates = datesSet.swaps;
swapsRates = ratesSet.swaps;
swapsRates = 0.5 * (swapsRates(:,1) + swapsRates(:,2));
% fix settlement date
t0 = datesSet.settlement;

% interpolate the first discount factor (swap 1y) and save it
delta = yearfrac(t0, swapsDates(1), EU_30_360);
discFact = intExtDF(discounts, dates, swapsDates(1));
% BPV represents BPV(0, T_i-1)
% initialize with the first Discount factor
BPV = delta * discFact;

% compute the discount factors for the swaps
for i=SwapStart:length(swapsDates)
    % compute delta for t_i-1, t_i
    delta = yearfrac(swapsDates(i-1), swapsDates(i), EU_30_360);
    % compute new discount factor for t_i
    discount = (1 - swapsRates(i) * BPV) / (1 + swapsRates(i) * delta);
    % update the dates, discounts and zero rates
    dates = [dates; swapsDates(i)];
    discounts = [discounts; discount];
    % update the BPV
    BPV = BPV + delta * discount;
end

end
