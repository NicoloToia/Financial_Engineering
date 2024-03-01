function [discount]=intExtDF(discounts, dates, targetDate)
% intExtDF: Interpolates (linear) or extrapolates (flat) the zero rates curve for a given date,
%           if the discount is present for the given date it returns the
%           value without interpolating or extrapolating
%
% INPUT
% discounts  : discount factors curve for preceding values (Bootstrap) 
% dates      : dates of discount factors curve bootstrap
% targetDate : corrisponding date of the discount requested
%
% OUTPUT
% discount   : discont factor found by matching/interpolating/extrapolating
%              the zero rates curve


% Zero rates are interpolated with yearfrac Convention ACT/365
ACT_365 = 3;

% Perfect match, return the corresponding discount
idx = find(dates==targetDate);
if ~isempty(idx)
    discount = discounts(idx);
    return;
end

% Find dates that encapsulate the target date
prevIdx = find(dates<targetDate,1,'last');
nextIdx = find(dates>targetDate,1,'first');

% Find the discount by linear interpolation of zero rates curve
if ~isempty(prevIdx) && ~isempty(nextIdx)
    % Set dates
    prevDate = dates(prevIdx);
    nextDate = dates(nextIdx);
    % Compute the zero rates
    prevY = -log(discounts(prevIdx)) / yearfrac(dates(1),prevDate, ACT_365);
    nextY = -log(discounts(nextIdx)) / yearfrac(dates(1),nextDate, ACT_365);
    % Linear interpolation between prevY & nextY and compute discount
    dt = yearfrac(prevDate,nextDate, ACT_365);
    y = prevY + (nextY - prevY) / dt * yearfrac(prevDate,targetDate, ACT_365);
    discount = exp(-y * yearfrac(dates(1),targetDate, ACT_365));
    return;
end

% Target date is beyond the last date, flat extrapolation
if targetDate > dates(end)
    discount = discounts(end);
    return;
end

end