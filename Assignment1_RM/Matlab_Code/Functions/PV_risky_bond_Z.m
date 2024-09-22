function PV = PV_risky_bond_Z(z, cf_schedule, ZC_curve)
% Compute the present value of a risky bond using the zero curve
% and the Z-score
%
% INPUTS
% zScore: the Z-score of the bond
% couponSchedule: a matrix with the coupon dates and values
% ZC_curve: a matrix with the zero curve
%
% OUTPUT
%   PV : Dirty price for a given risky bond (from scalar Z-spread)

% save the coupon dates and values
couponDates = cf_schedule(:,1);
couponValues = cf_schedule(:,2);

% compute the DFs
zeroDates = ZC_curve(:,1);
zeroRates = ZC_curve(:,2);
zeroRates = interp1(zeroDates, zeroRates, couponDates);
DF = exp(-zeroRates .* couponDates);
% compute the defaultable bond price
B_hat = DF .* exp(-z * couponDates);

% compute price by discounting the cash flows using the defaultable bond
PV = couponValues' * B_hat;

end