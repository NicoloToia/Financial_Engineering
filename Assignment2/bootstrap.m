function [dates, discounts, zeroRates]=bootstrap(datesSet, ratesSet)
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
ACT_365 = 3;
EU_30_360 = 6;

deposDates = datesSet.depos(1:3); % first 3 deposit dates
futuresDates = datesSet.futures(1:7,:); % first 7 futures dates
swapsDates = datesSet.swaps(2:end); % all swap dates except the first one

deposRates = ratesSet.depos(1:3,:) % first 3 deposit rates
futuresRates = ratesSet.futures(1:7,:); % first 7 futures rates
swapsRates = ratesSet.swaps(:,:); % all swap rates except the first one

% take the mid of the ask and bid
t0 = datesSet.settlement;
deposRates = 0.5 * (deposRates(:,1) + deposRates(:,2));
futuresRates = 0.5 * (futuresRates(:,1) + futuresRates(:,2));
swapsRates = 0.5 * (swapsRates(:,1) + swapsRates(:,2));

% create a vector of dates and discount factors
dates = [t0];
discounts = [1];
zeroRates = [0];

% use the act/360 convention for deposit rates
depoDelta = yearfrac(t0, deposDates, ACT_360)
dates = [dates; deposDates];
discounts = [discounts; 1 ./ (1 + deposRates .* depoDelta)]
zeroRates = [zeroRates; -log(discounts(2:end)) ./ depoDelta];

% futures rates (use act/360)
for i = 1:7
    % compute the forward discount factor
    delta = yearfrac(futuresDates(i, 1), futuresDates(i, 2), ACT_360);
    forwardDiscount = 1 / (1 + futuresRates(i) * delta);
    % compute the discount factor at the settlement date
    settlementDF = interpDF(discounts,zeroRates, dates, futuresDates(i, 1));
    % compute the discount factor at the expiry date
    discounts = [discounts; forwardDiscount * settlementDF];
    % update the zero rate
    delta = yearfrac(t0, futuresDates(i, 2), ACT_360);
    zeroRates = [zeroRates; -log(discounts(end)) / delta];
    dates = [dates; futuresDates(i, 2)];
end

% swaps rates (use act/365)
delta = yearfrac(t0, swapsDates(1), EU_30_360);
discFact = interpDF(discounts, zeroRates, dates, swapsDates(1));
BPV = delta * discFact;
for i=2:length(swapsDates)
    % compute new delta
    delta = yearfrac(swapsDates(i-1), swapsDates(i), EU_30_360);
    % compute new discount factor
    discount = (1 - swapsRates(i) * BPV) / (1 + swapsRates(i) * delta);
    % update the zero rate
    dates = [dates; swapsDates(i)];
    discounts = [discounts; discount];
    % update the BPV
    BPV = BPV + delta * discount;
end

% update the zero rates
zeroRates = -log(discounts) ./ yearfrac(t0, dates, EU_30_360);

end