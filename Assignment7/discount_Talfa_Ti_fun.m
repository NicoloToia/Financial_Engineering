function discount_Talfa_Ti = discount_Talfa_Ti_fun(alfa,omega,discounts)

    discount_Talfa_Ti = discounts(alfa+1:omega) / discounts(alfa);

end