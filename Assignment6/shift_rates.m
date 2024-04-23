function [rates_Set_shift] = shift_rates(ratesSet, datesSet, bucketDate,bp)

    % shift rates by 1bp only in the bucket date
    rates_Set_shift = ratesSet;

    dates = [
    datesSet.settlement;
    datesSet.depos(1:DeposDate);
    datesSet.futures(1:FutureDate,2);
    datesSet.swaps(SwapDate:end)
];

    % extract position of the bucket date from the datesSet
    [~,idx] = ismember(bucketDate, dates);

    % shift the rates
    rates_Set_shift(idx) = ratesSet(idx) + bp;

end