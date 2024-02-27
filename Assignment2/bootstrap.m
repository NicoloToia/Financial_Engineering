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

    %   disount and dates for Depos
    [discounts, dates] = discount_Depos(datesSet, ratesSet, DeposDate);
    %   disount and dates for Futures
    [discounts dates] =Discount_futures(datesSet, ratesSet, discounts, dates, FutureDate);
    %   disount and dates for Swaps
    [discounts dates] =Discount_swaps(datesSet, ratesSet, discounts, dates , SwapStart);

end