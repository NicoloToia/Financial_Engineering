function Jamshidian(alfa,omega,strike_swaption, node_dates, a, sigma, dates,discounts)

ACT_365 = 3;

t0 = dates(1);
Talfa = node_dates(alfa);

j=1;

discount_alfa = discounts(alfa);

Coupon_Call_option = 0;

ZC_call_option = zeros(1,omega-alfa-1);

for i = alfa+1 : omega

    discount_Talfa_Ti = discount_Talfa_Ti_fun(alfa,omega,discounts);

    strike_i = find_strike(alfa,omega,strike_swaption, Talfa, node_dates, discount_Talfa_Ti, a, sigma, t0);

    delta = yearfrac(Talfa,node_dates(i),ACT_365);
    
    % for every Ti compute the ZC call option 
    ZC_call_option = ZC_call(strike_i,discount_Talfa_Ti',discount_alfa,Talfa,t0,sigma,a,node_dates(alfa+1:omega));

    % coupon call option is the sum of the ZC call option for every Ti
    Coupon_Call_option = strike_swaption * delta * ZC_call_option(j) + Coupon_Call_option;

    j=j+1;
end

end