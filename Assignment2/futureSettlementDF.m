function [discount]=futureSettlementDF(discounts, dates, targetDate)

% perfect match, return the corresponding discount
idx = find(dates==targetDate);
if ~isempty(idx)
    discount = discounts(idx);
    return;
end

% there are dates that encapsulate the target date, interpolate
prevIdx = find(dates<targetDate,1,'last');
nextIdx = find(dates>targetDate,1,'first');

if ~isempty(prevIdx) && ~isempty(nextIdx)
    % get the dates
    prevDate = dates(prevIdx);
    nextDate = dates(nextIdx);
    % get the zero rates
    prevY = -log(discounts(prevIdx)) / yearfrac(dates(1),prevDate, 2);
    nextY = -log(discounts(nextIdx)) / yearfrac(dates(1),nextDate, 2);
    dt = yearfrac(prevDate,nextDate, 2);
    y = prevY + (nextY - prevY) / dt * yearfrac(prevDate,targetDate, 2);
    discount = exp(-y * yearfrac(dates(1),targetDate, 2));
    return;
end

% target date is beyond the last date, extrapolate
if targetDate > dates(end)
    % get last two dates and zero rates
    prevDate = dates(end-1);
    nextDate = dates(end);
    prevY = -log(discounts(end-1)) / yearfrac(dates(1),prevDate, 2);
    nextY = -log(discounts(end)) / yearfrac(dates(1),nextDate, 2);
    dt = yearfrac(prevDate,nextDate, 2);
    y = prevY + (nextY - prevY) / dt * yearfrac(prevDate,targetDate, 2);
    discount = exp(-y * yearfrac(dates(1),targetDate, 2));
    return;
end

end
 

