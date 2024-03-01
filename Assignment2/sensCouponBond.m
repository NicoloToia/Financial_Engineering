function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
% sensCouponBond compute the Macaulay duration
%
% INPUT
% setDate              : settlement date
% couponPaymentDates   : fixed coupon date of the bond
% fixedRate            : given fixed rate used for calculation
% dates                : Dates of the discount factors
% discounts            : Discount factors computed in the bootstrap

% OUTPUT
% MacD   : variation in the value of an IB coupon bond for
%          the parallel shift of the corresponding zero rate curve of 1bp,
%          normalized by bond price

% Compute the coupons to match yield of the swap
% Yearfrac Convention
EU_30_360 = 6;

% delta betweens coupons dates
deltas = [
    yearfrac(setDate, couponPaymentDates(1), EU_30_360);
    yearfrac(couponPaymentDates(1:end-1), couponPaymentDates(2:end), EU_30_360);
];
% compute the coupons & add the face falue of 1 to the last coupon
c = ones(length(couponPaymentDates),1) * fixedRate .* deltas;
c(end) = 1 + c(end);

% Compute the discounts factors
discountCoupons = zeros(length(couponPaymentDates),1);
for i=1:length(couponPaymentDates)
    % intExtDF is a function that interpolates the discount factors if necessary
    discountCoupons(i) = intExtDF(discounts, dates, couponPaymentDates(i));
end

% compute Macaulay duration
MacD = sum(c .* discountCoupons .* yearfrac(setDate, couponPaymentDates, EU_30_360)) / ...
    sum(c .* discountCoupons);
end