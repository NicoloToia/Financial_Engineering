function [discounts] = discountDepos(datesSet, ratesSet, discounts, DeposDate)
% Bootstrap the discount factors for deposits
%
% INPUT
% datesSet  : dates from market data stored as follow
%              -> datesSet.settlement : Settlement Date
%              -> datesSet.depos      : Vector of end calculation dates for Depos
%              -> datesSet.futures    : Matrix with start/end dates for Futures
%              -> datesSet.swaps      : Vector of end calculation date for Swaps
% ratesSet  : market data for bid&ask rates stored as follow
%              -> ratesSet.depos      : Matix bid&ask Depos
%              -> ratesSet.futures    : Matix bid&ask Futures
%              -> ratesSet.swaps      : Matix bid&ask Swaps
% discounts : Inizialized vector for discounts of the length needed, 
%             set as first value one
% DeposDate : the function compute the discounts for the first n depos,
%             where n = DeposDate
%
% OUTPUT
% discounts : Discount Factors for depos

% Yearfrac Convention 
ACT_360 = 2;

% Retrieve rates and dates
deposDates = datesSet.depos(1:DeposDate);               % first 3 depos dates
deposRates = ratesSet.depos(1:DeposDate,:);             % first 3 depos rates
deposRates = 0.5 * (deposRates(:,1) + deposRates(:,2)); % mid-market depos rates

% Fix settlement date
t0 = datesSet.settlement;

% Compute the discount factors for depos
depoDeltas = yearfrac(t0, deposDates, ACT_360);
discounts(2:DeposDate+1) = 1 ./ (1 + deposRates .* depoDeltas);

end