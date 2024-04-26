function shifted_ratesSet = shift_rate(ratesSet, datesSet, target_date, shift)
% shift_rate: shift the rate at the target date by a given amount
%
% INPUT
% ratesSet      : market data for mid-market rates stored as follow
%                 -> ratesSet.depos      : Matix bid&ask Depos
%                 -> ratesSet.futures    : Matix bid&ask Futures
%                 -> ratesSet.swaps      : Matix bid&ask Swaps
% datesSet     : dates of the market data stored as follow
%                 -> datesSet.settlement : Settlement Date
%                 -> datesSet.depos      : Vector of end calculation dates for Depos
%                 -> datesSet.futures    : Matrix with start/end dates for Futures
%                 -> datesSet.swaps      : Vector of end calculation date for Swaps
% target_date   : date at which the rate should be shifted
% shift         : amount by which the rate should be shifted

% find the number of fields in the struct
num_fields = numel(fieldnames(datesSet));
rates_names = fieldnames(ratesSet);
dates_names = fieldnames(datesSet);

% search for the target date in the datesSet
target_idx = NaN;
target_field = NaN;

for i = 2:num_fields
    for j = 1:length(datesSet.(dates_names{i}))
        target_row = datesSet.(dates_names{i})(j, :);
        if ismember(target_date, target_row)
            target_idx = j;
            target_field = i-1; % ratesSet has one less field than datesSet
            break
        end
    end
end

if isnan(target_idx)
    % print error message with the target date
    error(['The target date ', datestr(target_date), ' is not present in the datesSet.']);
end

% shift the rate at the target date by the given amount
shifted_ratesSet = ratesSet;
shifted_ratesSet.(rates_names{target_field})(target_idx) = ...
    ratesSet.(rates_names{target_field})(target_idx) + shift;

end