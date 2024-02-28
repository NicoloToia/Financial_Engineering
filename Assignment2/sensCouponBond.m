 function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)

    price=0;
    num=0;
    EU_30_360 = 6;
    deltas = [
        yearfrac(setDate, couponPaymentDates(1), EU_30_360);
        yearfrac(couponPaymentDates(1:end-1), couponPaymentDates(2:end), EU_30_360);
    ];
    c = ones(length(couponPaymentDates),1) * fixedRate .* deltas;
    c(end) = c(end) + 1;

    for i=1:length(couponPaymentDates)
        num = num+c(i)*intExtDF(discounts, dates, couponPaymentDates(i)) * yearfrac(setDate, couponPaymentDates(i), EU_30_360);
        price = price+c(i)*intExtDF(discounts, dates, couponPaymentDates(i));
    end
    MacD=num/price;
 end