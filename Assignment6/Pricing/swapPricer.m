function S = swapPricer(NPV, swapDates, discounts, dates)
% SWAPPRICER computes the fixed rate of a payer swap given the NPV
%
% INPUTS
%   NPV: NPV of the swap
%   swapDates: dates of the swap (including the settlement date)
%   discounts: discounts
%   dates: dates of the discounts

% discounts on the swap dates
discounts_swap = intExtDF(discounts, dates, swapDates(2:end));

% libor payments
libor_leg = discounts_swap(1) - discounts_swap(end);

% BPV
ACT_365 = 3;
deltas = yearfrac(swapDates(1:end-1), swapDates(2:end), ACT_365);

% BPV
BPV = deltas' * discounts_swap;

% find the correct fixed rate
S = (libor_leg - NPV) / BPV;

end