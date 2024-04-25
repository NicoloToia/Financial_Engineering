function [vols_shift] = shift_vols(spotVols, datesSet, bucketDate, percentage)

    % shift vols by 1% only in the bucket date
    vols_shift = spotVols;
    DeposDate = 3; 
    FutureDate = 7;
    SwapDate = 2; 
    dates = [
    datesSet.settlement;
    datesSet.depos(1:DeposDate);
    datesSet.futures(1:FutureDate,2);
    datesSet.swaps(SwapDate:end)];

    dates = datetime(dates, 'ConvertFrom', 'datenum');

    % extract position of the bucket date from the datesSet
    position = ismember(bucketDate, dates);
    idx = find(position);
    % shift the vols
    vols_shift(idx) = spotVols(idx) + percentage;

end