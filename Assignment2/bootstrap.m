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

depoDates = datesSet.depos(1:3); % first 3 deposit dates
futuresDates = datesSet.futures(1:7); % first 7 futures dates
swapDates = datesSet.swaps(2:end); % all swap dates except the first one

depoRates = ratesSet.depos(1:3); % first 3 deposit rates
futuresRates = ratesSet.futures(1:7); % first 7 futures rates
swapRates = ratesSet.swaps(2:end); % all swap rates except the first one

% take the mid of the ask and bid
depoRates = 0.5 * (depoRates(:,1) + depoRates(:,2));
futuresRates = 0.5 * (futuresRates(:,1) + futuresRates(:,2));
swapRates = 0.5 * (swapRates(:,1) + swapRates(:,2));

% create a vector of dates and discount factors
dates = [depoDates; futuresDates(:, 2); swapDates];
discounts = zeros(size(dates));

% deposit rates (use act/360 convention)
depoDeltas = yearfrac(depoDates(1), depoDates(2), 2);
discounts(1:3) = 1 ./ (1 + depoRates .* depoDeltas);

% futures rates (use act/360 convention)
for i=4:10
    % compute the forward discount factor
    delta = yearfrac(futuresDates(i, 1), futuresDates(i, 2), 2);
    forwardDiscount = 1 / (1 + futuresRates(i) * delta);
    % find the index of the previous and next dates
    prevIdx = find(dates <= futuresDates(i, 1), 1, 'last');
    nextIdx = find(dates >= futuresDates(i, 2), 1, 'first');

    % compute the discount factor for the settlement date
    settlementDF = interpDF(futuresDates(i, 1), dates(prevIdx), dates(nextIdx), discounts(prevIdx), discounts(nextIdx));

    % compute the discount factor for the expiry date
    discounts(i) = forwardDiscount * settlementDF;

end

end

