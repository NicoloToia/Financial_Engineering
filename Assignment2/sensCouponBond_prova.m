 function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
    c=ones(length(couponPaymentDates),1)*fixedRate;
    % find the appropriate discount factors
    idx = ismember(dates, couponPaymentDates);
    discountsPaymentDates = discounts(idx);
    % transform the dates in years
    t = yearfrac(setDate, couponPaymentDates, 6);

    MacD = sum(c.*t.*discountsPaymentDates)/sum(c.*discountsPaymentDates);
 end