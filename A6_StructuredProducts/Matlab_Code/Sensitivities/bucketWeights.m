function bucket_sens = bucketWeights(settlement_date, sensitivities, sens_dates, start, center, end_, first, last)
% compute the buckets for the center date
%
% INPUT
%   sensitivities : sensitivities
%   dates : dates of the sensitivities
%   start : start date of the sensitivities (in years)
%   center : center date of the sensitivities (in years)
%   end : end date of the sensitivities (in years)
%   first : first bucket end date (in years)
%   last : last bucket end date (in years)

% first bucket 0-2y

% compute the center and last date
begin_date = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(start);
center_date = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(center);
end_date = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(end_);
% move to business days
if ~isbusday(begin_date, eurCalendar())
    begin_date = busdate(begin_date, 'modifiedfollow', eurCalendar());
end
if ~isbusday(center_date, eurCalendar())
    center_date = busdate(center_date, 'modifiedfollow', eurCalendar());
end
if ~isbusday(end_date, eurCalendar())
    end_date = busdate(end_date, 'modifiedfollow', eurCalendar());
end

% convert to datenum
begin_date = datenum(begin_date);
center_date = datenum(center_date);
end_date = datenum(end_date);

% find the relevant dates
first_half_dates = sens_dates(sens_dates >= begin_date & sens_dates < center_date);
first_half_sens = sensitivities(sens_dates >= begin_date & sens_dates < center_date);

center_sens = sensitivities(sens_dates == center_date);
center_date = sens_dates(sens_dates == center_date);

second_half_sens = sensitivities(sens_dates > center_date & sens_dates <= end_date);
second_half_dates = sens_dates(sens_dates > center_date & sens_dates <= end_date);

% compute the weights
weights = [
    (first_half_dates - begin_date) / (center_date - begin_date);
    1;
    (end_date - second_half_dates) / (end_date - center_date);
];

% adjust for the first and last bucket
if first
    % adjust the weights (all ones up to the center date)
    weights(1:length(first_half_dates)) = 1;
end

if last
    % adjust the weights (all ones after the center date)
    weights(length(first_half_dates)+1:end) = 1;
end

% compute the bucket sensitivity
bucket_sens = weights' * [first_half_sens; center_sens; second_half_sens];

end