 function MacD = sensCouponBond_prova(setDate, couponPaymentDates, fixedRate, dates, discounts)

    MacD=0;
    price=0;
    num=0;
    c=ones(length(couponPaymentDates),1)*fixedRate;
    c(length(couponPaymentDates))=fixedRate+1;

    for i=1:length(couponPaymentDates)
        num = num+c(i)*(dates(i)-setDate)*futureSettlementDF(discounts, dates, couponPaymentDates(i));
        price = price+c(i)*futureSettlementDF(discounts, dates, couponPaymentDates(i));
    end
    MacD=num/price;
 end