function NPV = swapNPV(swapRate, ttm, discounts, dates)
% SWAPNPV computes the NPV of a payer swap with the given fixed rate
%
% INPUTS
%   swapRate: fixed rate of the swap
%   ttm: time to maturity of the swap (in years)
%   discounts: discounts
%   dates: dates of the discounts

% compute the dates
swapDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(0:ttm)';

% move to business days if needed
swapDates(~isbusday(swapDates, eurCalendar())) = ...
    busdate(swapDates(~isbusday(swapDates, eurCalendar())), 'modifiedfollow', eurCalendar());
swapDates = datenum(swapDates);

% discounts on the swap dates
discounts_swap = intExtDF(discounts, dates, swapDates(2:end));

% libor payments
libor_leg = 1 - discounts_swap(end);

% BPV
EU_30_360 = 6;
deltas = yearfrac(swapDates(1:end-1), swapDates(2:end), EU_30_360);

% BPV
BPV = deltas' * discounts_swap;

% fixed leg
fixed_leg = swapRate * BPV;

% NPV
NPV = libor_leg - fixed_leg;

end