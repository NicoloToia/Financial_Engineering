function [discount]=interpDF(discounts,zeroRates, dates, targetDate)

% find the date that matches the target date
idx = find(dates==targetDate);
if ~isempty(idx)
    discount = discounts(idx);
    return;
end

% if there are dates that encapsulate the target date, interpolate
prevIdx = find(dates<targetDate,1,'last');
nextIdx = find(dates>targetDate,1,'first');

if ~isempty(prevIdx) && ~isempty(nextIdx)
    prevDate = dates(prevIdx);
    nextDate = dates(nextIdx);
    prevY = zeroRates(prevIdx);
    nextY = zeroRates(nextIdx);
    delta = yearfrac(prevDate,nextDate, 2);
    y = prevY + (nextY - prevY) / delta * yearfrac(prevDate,targetDate, 2);
    discount = exp(-y * yearfrac(dates(1),targetDate, 2));
    return;
end

% if the date is beyond the last date, extrapolate
if targetDate > dates(end)
    prevDate = dates(end-1);
    nextDate = dates(end);
    prevY = zeroRates(end-1);
    nextY = zeroRates(end);
    delta = yearfrac(prevDate,nextDate, 2);
    y = prevY + (nextY - prevY) / delta * yearfrac(prevDate,targetDate, 2);
    discount = exp(-y * yearfrac(dates(1),targetDate, 2));
    return;
end

end
 

