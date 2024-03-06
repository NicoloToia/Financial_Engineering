function PV = PV_risky_bond_Z(zScore, couponSchedule, ZC_curve)

% save the coupon dates and values
couponDates = couponSchedule(:,1);
couponValues = couponSchedule(:,2);

% compute the DFs
zeroDates = ZC_curve(:,1);
zeroRates = ZC_curve(:,2);
zeroRates = interp1(zeroDates, zeroRates, couponDates);
DF = exp(-zeroRates .* couponDates);

B_hat = DF .* exp(-zScore * couponDates);

PV = couponValues' * B_hat;

end