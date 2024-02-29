function [DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, ...
    discounts,discounts_DV01)

% 1 bp is 0.01% or 0.0001
bp = 0.0001;

% find the discounts factors for the fixed leg
discountsFixedLeg = zeros(length(fixedLegPaymentDates), 1);
discountsFixedLeg_DV01 = zeros(length(fixedLegPaymentDates), 1);
for i = 1:length(fixedLegPaymentDates)
    % intExtDF is a function that interpolates the discount factors if necessary
    discountsFixedLeg(i) = intExtDF(discounts, dates, fixedLegPaymentDates(i));
    discountsFixedLeg_DV01(i) = intExtDF(discounts_DV01, dates, fixedLegPaymentDates(i));
end

% compute the BPV
deltas = [
    yearfrac(setDate, fixedLegPaymentDates(1));
    yearfrac(fixedLegPaymentDates(1:end-1), fixedLegPaymentDates(2:end));
];
BPV = deltas' * discountsFixedLeg * bp;

% original NPV
NPV_0 = 1 - discountsFixedLeg(end) - fixedRate * deltas' * discountsFixedLeg;

% compute the DV01
NPV_shift = 1 - discountsFixedLeg_DV01(end) - fixedRate * deltas' * discountsFixedLeg_DV01;
DV01 = NPV_shift - NPV_0;

% parallel shift of zero rates
zeroRatesFixedLeg = zeroRates([setDate;fixedLegPaymentDates], [1;discountsFixedLeg])/100 + bp;
zeroRatesFixedLeg = zeroRatesFixedLeg(2:end); % remove the zero rate at the set date
discountsFixedLeg_z = exp(-zeroRatesFixedLeg .* yearfrac(setDate, fixedLegPaymentDates));

% compute the NPV with the shifted zero rates
NPV_z = 1 - discountsFixedLeg_z(end) - fixedRate * deltas' * discountsFixedLeg_z;
DV01_z = NPV_z - NPV_0;

end