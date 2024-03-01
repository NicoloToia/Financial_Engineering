function [dates, discounts]=bootstrap(datesSet, ratesSet)
% Bootstrap the discount factors curve from the input dates and rates
%
% INPUT
% datesSet: dates from market data stored as follow
%            -> datesSet.settlement : Settlement Date
%            -> datesSet.depos      : Vector of end calculation dates for Depos
%            -> datesSet.futures    : Matrix with start/end dates for Futures
%            -> datesSet.swaps      : Vector of end calculation date for Swaps
% ratesSet: market data for bid&ask rates stored as follow
%            -> ratesSet.depos      : Matix bid&ask Depos
%            -> ratesSet.futures    : Matix bid&ask Futures
%            -> ratesSet.swaps      : Matix bid&ask Swaps
%
% OUTPUT
% dates    : Dates of discount factors curve bootstrap
%            -> recall that the bootstrap is considering the most liquid
%            products, hence as default & "according" with the market the vector has the 
%            following structure:
%            [settelment date; first 3 depos expiries; first 7 futures expiries; from the 2nd swap expiry to end]
% discounts: Discount Factors of the discount curve bootstrap at corresponding dates

% Parameters fixed by trader, based on liquidity
DeposDate  = 3;     %   consider until this Depos
FutureDate = 7;     %   consider until this Future
SwapDate  = 2;     %   consider from this Swap

% Dates vector initialization thanks to parameters above
dates = [
    datesSet.settlement;
    datesSet.depos(1:DeposDate);
    datesSet.futures(1:FutureDate,2);
    datesSet.swaps(SwapDate:end)
];

% Initialize discounts vector
discounts = ones(length(dates),1);

% indices
futureStart = 1 + DeposDate;
futureEnd = futureStart + FutureDate;
swapStart = futureEnd + 1;

% discount for Depos
[discounts] = discountDepos(datesSet, ratesSet, discounts, DeposDate);
% discount for Futures
[discounts] = discountFutures(datesSet, ratesSet, discounts, dates, FutureDate, futureStart);
% discount for Swaps
[discounts] = discountSwaps(datesSet, ratesSet, discounts, dates , SwapDate, swapStart);

end
