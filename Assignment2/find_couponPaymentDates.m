function [couponDates] = find_couponPaymentDates(datesSet, setDate)
    %   coupon payments every 3 months for 6 years
    couponDates=zeros(24,1);
    for i=1:24
        couponDates(i)=setDate+90;
    end
end