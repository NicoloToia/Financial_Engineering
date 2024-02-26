function [dates, discounts]=bootstrap(datesSet, ratesSet)
% Bootstrap the discount factors curve from the input dates and rates
%
% INPUTS:
% datesSet: a vector of dates
% ratesSet: a vector of rates
%
% OUTPUTS:
% dates: a vector of dates
% discounts: a vector of discount factors

% dates conventions
ACT_360 = 2;
EU_30_360 = 6;

% settlement date
t0 = datesSet.settlement;
% first 3 deposits
deposDates = datesSet.depos(1:3);
deposRates = ratesSet.depos(1:3,:)
% first 7 futures
futuresDates = datesSet.futures(1:7,:);
futuresRates = ratesSet.futures(1:7,:);
% all swaps (first will be ignored later)
swapsDates = datesSet.swaps;
swapsRates = ratesSet.swaps;

% take the mid of the ask and bid
deposRates = 0.5 * (deposRates(:,1) + deposRates(:,2));
futuresRates = 0.5 * (futuresRates(:,1) + futuresRates(:,2));
swapsRates = 0.5 * (swapsRates(:,1) + swapsRates(:,2));

% initialize dates and discounts
dates = [t0];
discounts = [1];

% Deposit rates
% use act/360
depoDeltas = yearfrac(t0, deposDates, ACT_360)
dates = [dates; deposDates];
discounts = [discounts; 1 ./ (1 + deposRates .* depoDeltas)]

% futures rates (use act/360)
for i = 1:7
    % compute the forward discount factor
    deltaForward = yearfrac(futuresDates(i, 1), futuresDates(i, 2), ACT_360);
    forwardDiscount = 1 / (1 + futuresRates(i) * deltaForward);
    % compute the discount factor at the settlement date
    settlementDF = futureSettlementDF(discounts, dates, futuresDates(i, 1));
    % TODO: add the settlement discount and date to the vectors
    % compute the discount factor at the expiry date
    discounts = [discounts; forwardDiscount * settlementDF];
    % update the zero rate
    delta = yearfrac(t0, futuresDates(i, 2), ACT_360);
    dates = [dates; futuresDates(i, 2)];
end

% swaps rates (use EU 30/360)
delta = yearfrac(t0, swapsDates(1), EU_30_360);
% interpolate the first discount factor
discFact = futureSettlementDF(discounts, dates, swapsDates(1));
% BPV represents BPV(0, T_i-1)
% initialize with the first Discount factor
BPV = delta * discFact;

for i=2:length(swapsDates)
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