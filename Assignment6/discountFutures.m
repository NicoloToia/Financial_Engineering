function [discounts] = discountFutures(datesSet, ratesSet, discounts, dates, FutureDate, futureStart)
% Bootstrap the discount factors for futures
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
% FutureDate : the function compute the discounts for the first n futures,
%              where n = FutureDate
%
% OUTPUT
% discounts: Discount Factors for depos
    
% Yearfrac Convention 
ACT_360 = 2;

% Retrieve rates and dates
futuresDates = datesSet.futures(1:FutureDate,:);              % first 7 futures dates
futuresRates = ratesSet.futures(1:FutureDate);              % first 7 futures rates

% Compute the discount factors for the futures
for i = 1:FutureDate
    % compute the forward discount factor
    deltaForward = yearfrac(futuresDates(i, 1), futuresDates(i, 2), ACT_360);
    forwardDiscount = 1 / (1 + futuresRates(i) * deltaForward);
    % find the discount factor at the settlement date by intExDF functions
    settlementDF = intExtDF(discounts(1:futureStart + i -1), dates(1:futureStart + i - 1), futuresDates(i, 1));
    % compute the discount factor at the expiry date
    discounts(futureStart + i) = forwardDiscount * settlementDF;
end

end