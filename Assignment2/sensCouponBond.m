 function MacD = sensCouponBond(setDate, couponPaymentDates, fixedRate, dates, discounts)
    MacD=0;
    c=ones(length(couponPaymentDates),1)*fixedRate;
    c(length(couponPaymentDates))=fixedRate+1;
    for i=1:length(couponPaymentDates)
        if dates(i)=couponPaymentDates(i)
            MacD=MacD+c(i)*(dates(i)-setDate)*discounts(i)/(c(i)*discounts(i))
        else
            MacD=MacD+c(i)*(dates(i)-setDate)*futureSettlementDF(discounts, dates, couponPaymentDates(i))/(c(i)*futureSettlementDF(discounts, dates, couponPaymentDates(i)))
        end
    end
 end