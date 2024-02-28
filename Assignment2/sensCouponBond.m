function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)

% compute the coupons to match yield of the swap
EU_30_360 = 6;
deltas = [
    yearfrac(setDate, couponPaymentDates(1), EU_30_360);
    yearfrac(couponPaymentDates(1:end-1), couponPaymentDates(2:end), EU_30_360);
];
c = ones(length(couponPaymentDates),1) * fixedRate .* deltas;
c(end) = 1 + c(end);

% compute the discounts factors
discountCoupons = zeros(length(couponPaymentDates),1);
for i=1:length(couponPaymentDates)
    discountCoupons(i) = intExtDF(discounts, dates, couponPaymentDates(i));
end

MacD = sum(c .* discountCoupons .* yearfrac(setDate, couponPaymentDates, EU_30_360)) / ...
    sum(c .* discountCoupons);
end