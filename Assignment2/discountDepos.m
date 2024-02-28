function [discounts, dates] = discountDepos(datesSet, ratesSet, DeposDate)

% date notation 
ACT_360 = 2;

% retrieve rates and dates
deposDates = datesSet.depos(1:DeposDate); % first 3 deposit dates
deposRates = ratesSet.depos(1:DeposDate,:); % first 3 deposit rates
deposRates = 0.5 * (deposRates(:,1) + deposRates(:,2));

% settlement date and initialize output
t0 = datesSet.settlement;
dates = [t0];
discounts = [1];
% compute the discount factors
depoDeltas = yearfrac(t0, deposDates, ACT_360);
dates = [dates; deposDates];
discounts = [discounts; 1 ./ (1 + deposRates .* depoDeltas)];

end