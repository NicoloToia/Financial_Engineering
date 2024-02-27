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

    % interpolate the remaining dates

    % depos
    for i=DeposDate+1:length(datesSet.depos)
        discount = futureSettlementDF(discounts, dates, datesSet.depos(i));
        dates = [dates;datesSet.depos(i)];
        discounts = [discounts; discount];
    end

    % futures
    for i=FutureDate+1:length(datesSet.futures)
        discount = futureSettlementDF(discounts, dates, datesSet.futures(i,2));
        dates = [dates;datesSet.futures(i,2)];
        discounts = [discounts; discount];
    end

    % swaps
    % interpolation uses ACT/360, fine since we are interpolating between
    % future dates
    discount = futureSettlementDF(discounts, dates, datesSet.swaps(1));
    dates = [dates;datesSet.swaps(1)];
    discounts = [discounts; discount];

    % sort by dates
    [dates, idx] = sort(dates);
    discounts = discounts(idx);

end