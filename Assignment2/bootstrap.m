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

deposDates = datesSet.depos(1:3); % first 3 deposit dates
futuresDates = datesSet.futures(1:7,:); % first 7 futures dates
swapsDates = datesSet.swaps(2:end); % all swap dates except the first one

deposRates = ratesSet.depos(1:3); % first 3 deposit rates
futuresRates = ratesSet.futures(1:7,:); % first 7 futures rates
swapsRates = ratesSet.swaps(2:end); % all swap rates except the first one

% take the mid of the ask and bid
t0 = datesSet.settlement;
deposRates = 0.5 * (deposRates(:,1) + deposRates(:,2));
futuressRates = 0.5 * (futuresRates(:,1) + futuresRates(:,2));
swapsRates = 0.5 * (swapsRates(:,1) + swapsRates(:,2));

% create a vector of dates and discount factors
dates = [t0];
discounts = [1];
zeroRates = [0];

% use the act/360 convention for deposit rates
depoDelta = yearfrac(t0, deposDates, 2);
dates = [dates, deposDates];
discounts = [discounts, 1 ./ (1 + deposRates .* depoDelta)];
zeroRates = [zeroRates, -log(discounts(2:end)) ./ depoDelta];

% futures rates (use act/360)
for i = 1:7
    % compute the forward discount factor
    delta = yearfrac(futuresDates(i, 1), futuresDates(i, 2), 2);
    forwardDiscount = 1 / (1 + futuresRates(i) * delta);
    % compute the discount factor at the settlement date
    settlementDF = interpDF(discounts,zeroRates, dates, futuresDates(i, 1));
    % compute the discount factor at the expiry date
    discounts(i+4) = forwardDiscount * settlementDF;
    % update the zero rate
    zeroRates(i+4) = -log(discounts(i+4)) / delta;
end

% % futures rates (use act/365 convention)
% for i=1:7
%     % compute the forward discount factor
%     delta = yearfrac(futuresDates(i, 1), futuresDates(i, 2), 2);
%     forwardDiscount = 1 / (1 + futuresRates(i) * delta);
%     % compute the discount factor for the settlement date
%     settlementDF = interpDF(dates, discounts, futuresDates(i, 2));
%     % compute the discount factor for the expiry date
%     discounts(i+length(depoRates)) = forwardDiscount * settlementDF;
% end

end

