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

%   parameters fixed by trader
DeposDate=3;    %   consider untill this Depo
FutureDate=7;    %   consider untill this Future
SwapStart=2;    %   consider ufrom this Swap

% discount and dates for Depos
[discounts dates] = discountDepos(datesSet, ratesSet, DeposDate);
% discount and dates for Futures
[discounts dates] = discountFutures(datesSet, ratesSet, discounts, dates, FutureDate);
% discount and dates for Swaps
[discounts dates] = discountSwaps(datesSet, ratesSet, discounts, dates , SwapStart);

end