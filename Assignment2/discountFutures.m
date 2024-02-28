function [discounts dates] = discountFutures(datesSet, ratesSet, discounts, dates, FutureDate)
    
% date notation 
ACT_360 = 2;
%   retrive data 
futuresDates = datesSet.futures(1:FutureDate,:); % first 7 futures dates
futuresRates = ratesSet.futures(1:FutureDate,:); % first 7 futures rates
futuresRates = 0.5 * (futuresRates(:,1) + futuresRates(:,2));

% compute rates from futures
for i = 1:FutureDate
    % compute the forward discount factor
    deltaForward = yearfrac(futuresDates(i, 1), futuresDates(i, 2), ACT_360);
    forwardDiscount = 1 / (1 + futuresRates(i) * deltaForward);
    % compute the discount factor at the settlement date
    settlementDF = futureSettlementDF(discounts, dates, futuresDates(i, 1));
    % compute the discount factor at the expiry date
    discounts = [discounts; forwardDiscount * settlementDF];
    % update the zero rate
    dates = [dates; futuresDates(i, 2)];
end

end