function S = swapPricer(NPV, ttm, discounts, dates)
% SWAPPRICER computes the fixed rate of a payer swap given the NPV
%
% INPUTS
%   NPV: NPV of the swap
%   ttm: time to maturity of the swap (in years)
%   discounts: discounts
%   dates: dates of the discounts

% compute the dates
swapDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(0:ttm)';
swapDates(~isbusday(swapDates, eurCalendar())) = ...
    busdate(swapDates(~isbusday(swapDates, eurCalendar())), 'modifiedfollow', eurCalendar());
swapDates
swapDates = datenum(swapDates);

% discounts on the swap dates
discounts_swap = intExtDF(discounts, dates, swapDates(2:end));

% libor payments
libor_leg = 1 - discounts_swap(end);

% BPV
ACT_365 = 3;
deltas = yearfrac(swapDates(1:end-1), swapDates(2:end), ACT_365)

% BPV
BPV = deltas' * discounts_swap;

% find the correct fixed rate
S = (libor_leg - NPV) / BPV;

end