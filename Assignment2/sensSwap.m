function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, ...
    discounts,discounts_DV01)
% sensSwap return the main sensitivities for a swap
%
% INPUT
% setDate              : settlement date
% fixedLegPaymentDates : fixed legs date of the swap (when S is paid or recieved)
% fixedRate            : given fixed rate used for calculation
% dates                : Dates of the discount factors
% discounts            : Discount factors computed in the bootstrap
% discounts_DV01       : Discount factors computed by the bootstrap of the shifted market data
%
% OUTPUT
% DV01   : Variation in the NPV for a derivative's portfolio (e.g. a signle
%          swap) due to the increase of rates used for the bootstrap of the
%          original curve for 1 bp in parallel shift
% BPV    : Net Present Value of 1 bp       
% DV01_z : variation in the NPV for a derivative's portfolio due to the
%          increase of zero rates for 1 bp in parallel shift; i.e. it is the
%          derivative's sensitivity for a (parallel shift)movement in the zero rate
%          curve    

% 1 bp is 0.01% or 0.0001
bp = 0.0001;

% find the discounts factors for the fixed leg
discountsFixedLeg = intExtDF(discounts, dates, fixedLegPaymentDates);
discountsFixedLeg_DV01 = intExtDF(discounts_DV01, dates, fixedLegPaymentDates);

% Compute the BPV
deltas = [
    yearfrac(setDate, fixedLegPaymentDates(1));
    yearfrac(fixedLegPaymentDates(1:end-1), fixedLegPaymentDates(2:end));
];
BPV = deltas' * discountsFixedLeg * bp;

% Original NPV (before the shift)
NPV_0 = 1 - discountsFixedLeg(end) - fixedRate * deltas' * discountsFixedLeg;
% Shifted NPV (after the shift of 1 bp)
NPV_shift = 1 - discountsFixedLeg_DV01(end) - fixedRate * deltas' * discountsFixedLeg_DV01;
% Compute the DV01
DV01 = abs(NPV_shift - NPV_0);

% Parallel shift of zero rates curve
zeroRatesFixedLeg = zeroRates([setDate;fixedLegPaymentDates], [1;discountsFixedLeg])/100 + bp;
zeroRatesFixedLeg = zeroRatesFixedLeg(2:end); % remove the zero rate at the settlemnet date
% Compute discounts
discountsFixedLeg_z = exp(-zeroRatesFixedLeg .* yearfrac(setDate, fixedLegPaymentDates));

% Compute the NPV of the shifted zero rates curve
NPV_z = 1 - discountsFixedLeg_z(end) - fixedRate * deltas' * discountsFixedLeg_z;
% Compute the DV01_z
DV01_z = abs(NPV_z - NPV_0);

end
