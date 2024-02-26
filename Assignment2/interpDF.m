function [discount]=interpDF(discounts,zeroRates, dates, targetDate)

% find the date that matches the target date
idx = find(dates==targetDate);
if ~isempty(idx)
    discount = discounts(idx);
    return;
end

% if there are dates that encapsulate the target date, interpolate










% find the index of the dates that encapsulate the target date
prevIdx = find(dates<=targetDate,1,'last');
nextIdx = find(dates>targetDate,1,'first');
prevDate = dates(prevIdx);
nextDate = dates(nextIdx);

% target dates coincides
if targetDate == prevDate
    y = zeroRates(prevIdx);
else
    prevRate = zeroRates(prevIdx)
    nextRate = zeroRates(nextIdx)
    y = prevRate + (nextRate - prevRate) / (nextDate - prevDate) * (targetDate - prevDate);
end

delta = yearfrac(dates(1),targetDate, 2);

discount = exp(-y*delta);

end
 

