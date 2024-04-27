function NPV = swapNPV(swapRate, swapDates, discounts, dates)
% SWAPNPV computes the NPV of a payer swap with the given fixed rate
%
% INPUTS
%   swapRate: fixed rate of the swap
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

% fixed leg
fixed_leg = swapRate * BPV;

% NPV
NPV = libor_leg - fixed_leg;

end