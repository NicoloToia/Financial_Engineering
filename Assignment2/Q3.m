% Price a 6y “I.B. coupon bond” issued on the 15th of Feb ’2008 with coupon rate equal to the 
% corresponding mid-market 7y swap rate. Assume for the coupons a 30/360 European day count. 
% Hint: Shortcuts are appreciated.

EU_Act = 6;
date_sw7 = find(dates == datenum('02/19/15'));

S = 0.5*(ratesSet.swaps(date_sw7-10,1) + ratesSet.swaps(date_sw7-10,2));

delta = yearfrac(dates(date_sw7-1),dates(date_sw7),EU_Act);

price_IBbond = (1 - discounts(date_sw7)) - S * delta * discounts(date_sw7) + discounts(date_sw7 - 1)

priceB = (1 + delta*S) * discounts(date_sw7) + discounts(date_sw7 - 1)

%%%%%%%%%%%%%

% find the discounts factors for the fixed leg
fixedLegPaymentDates = datesSet.swaps(1:6);

discountsFixedLeg = zeros(length(fixedLegPaymentDates), 1);
discountsFixedLeg_DV01 = zeros(length(fixedLegPaymentDates), 1);
for i = 1:length(fixedLegPaymentDates)
    % intExtDF is a function that interpolates the discount factors if necessary
    discountsFixedLeg(i) = intExtDF(discounts, dates, fixedLegPaymentDates(i));
    discountsFixedLeg_DV01(i) = intExtDF(discounts_DV01, dates, fixedLegPaymentDates(i));
end


couponPaymentDates = fixedLegPaymentDates;
EU_30_360 = 6;
deltas = [
    yearfrac(setDate, couponPaymentDates(1), EU_30_360);
    yearfrac(couponPaymentDates(1:end-1), couponPaymentDates(2:end), EU_30_360);
];
c = ones(length(couponPaymentDates),1) * S .* deltas;
c(end) = 1 + c(end);

% compute the discounts factors
discountCoupons = zeros(length(couponPaymentDates),1);
for i=1:length(couponPaymentDates)
    discountCoupons(i) = intExtDF(discounts, dates, couponPaymentDates(i));
end

priceB_classic = sum(c .* discountCoupons)
