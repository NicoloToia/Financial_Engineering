function z = zscore(ZC_curve, couponSchedule, dirtyPrice)

% save coupon dates and values
couponDates = couponSchedule(:,1);
couponValues = couponSchedule(:,2);

% compute relevant discount factors
zeroDates = ZC_curve(:,1);
zeroRates = ZC_curve(:,2);
zeroRates = interp1(zeroDates, zeroRates, couponDates);
DF = exp(-zeroRates .* couponDates);

P = @(z) couponValues' * (DF .* exp(-z * couponDates)) - dirtyPrice;

z = fzero(P, 0);

end 