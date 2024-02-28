function [discount]=intExtDF(discounts, dates, targetDate)
% intExtDF: Interpolates or extrapolates the discount factor for a given date
%
% discount = intExtDF(discounts, dates, targetDate) returns the discount
% factor for the given target date. The discount factors are given in the
% vector discounts and the corresponding dates are given in the vector
% dates.

% zero rates are interpolated with ACT/365
ACT_365 = 3;

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
    prevY = -log(discounts(prevIdx)) / yearfrac(dates(1),prevDate, ACT_365);
    nextY = -log(discounts(nextIdx)) / yearfrac(dates(1),nextDate, ACT_365);
    dt = yearfrac(prevDate,nextDate, ACT_365);
    y = prevY + (nextY - prevY) / dt * yearfrac(prevDate,targetDate, ACT_365);
    discount = exp(-y * yearfrac(dates(1),targetDate, ACT_365));
    return;
end

% target date is beyond the last date, extrapolate with flat
if targetDate > dates(end)
    discount = discounts(end);
    return;
end

end
 

