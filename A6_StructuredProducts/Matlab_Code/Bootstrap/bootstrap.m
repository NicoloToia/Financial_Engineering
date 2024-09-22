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

% all swap dates
swapDates = datetime(datesSet.settlement, 'ConvertFrom', 'datenum') + calyears(1:50)';
swapDates(~isbusday(swapDates, eurCalendar())) = busdate(swapDates(~isbusday(swapDates, eurCalendar())), 'modifiedfollow', eurCalendar());
swapDates = datenum(swapDates);

% Interpolate the mid market rates for swaps from the market data

% actually quoted swap dates
datesSet.swaps = [swapDates(1:12); swapDates(15:5:30); swapDates(40:10:50)];

% Interpolation of the mid market rates for swaps using spline
ACT_365 = 3;
delta_swaps_set = yearfrac(datesSet.settlement, datesSet.swaps, ACT_365);
delta_swaps = yearfrac(datesSet.settlement, swapDates, ACT_365);
ratesSet.swaps = interp1(delta_swaps_set, ratesSet.swaps, delta_swaps, 'spline');
datesSet.swaps = swapDates;

% Parameters fixed by trader, based on liquidity
DeposDate  = 4;     %   consider until this Depos
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
