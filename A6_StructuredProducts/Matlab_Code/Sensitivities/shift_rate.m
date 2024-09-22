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

% search for the target date in the datesSet
target_idx = NaN;

% search the depo dates
for i = 1:4
    if datesSet.depos(i) == target_date
        target_idx = i;
        break
    end
end

% target date found in the depos
if ~isnan(target_idx)
    shifted_ratesSet = ratesSet;
    shifted_ratesSet.depos(target_idx) = ratesSet.depos(target_idx) + shift;
    return
end

% search the futures dates
for i = 1:7
    if datesSet.futures(i, 2) == target_date
        target_idx = i;
        break
    end
end

% target date found in the futures
if ~isnan(target_idx)
    shifted_ratesSet = ratesSet;
    shifted_ratesSet.futures(target_idx) = ratesSet.futures(target_idx) + shift;
    return
else

% search the swaps dates
for i = 2:length(datesSet.swaps)
    if datesSet.swaps(i) == target_date
        target_idx = i;
        break
    end
end

% target date found in the swaps
if ~isnan(target_idx)
    shifted_ratesSet = ratesSet;
    shifted_ratesSet.swaps(target_idx) = ratesSet.swaps(target_idx) + shift;
    return
end

if isnan(target_idx)
    % print error message with the target date
    error(['The target date ', datestr(target_date), ' is not present in the datesSet.']);
end

end