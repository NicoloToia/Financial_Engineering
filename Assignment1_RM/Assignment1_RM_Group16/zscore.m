function [z] = zscore(ZC_curve, couponSchedule, dirtyPrice)
% Compute Z_scores
%
% INPUT
% zeroCurve : ZC bond data [yearfrac ; rates]
% couponSchedule: cash flows bond [yearfrac; cash flow]
% dirtyPrice = market dirty price
%
% OUTPUT
% z = z_score

% save coupon dates and values
couponDates = couponSchedule(:,1);
couponValues = couponSchedule(:,2);

% compute relevant discount factors
zeroDates = ZC_curve(:,1);
zeroRates = ZC_curve(:,2);
zeroRates = interp1(zeroDates, zeroRates, couponDates);
DF = exp(-zeroRates .* couponDates);

% Inizializate function in z
P = @(z) couponValues' * (DF .* exp(-z * couponDates)) - dirtyPrice;

% Find zero
z = fzero(P, 0);

end 