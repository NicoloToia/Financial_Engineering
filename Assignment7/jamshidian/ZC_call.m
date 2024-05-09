function ZC_call_price = ZC_call(strike_i,discount_Talfa_Ti,discount_alfa,Talfa,t0,sigma,a,dates_years)

ACT_365 = 3;

V = Volatility(Talfa,t0,dates_years,sigma,a)';

dt = yearfrac(t0,dates_years,ACT_365)';


% define d1 and d2 for every T_i given
d1 = log(discount_Talfa_Ti ./ strike_i)./(V.* sqrt(dt(1))) + 0.5 .* V .* sqrt(dt(1));
d2 = log(discount_Talfa_Ti ./ strike_i)./(V.* sqrt(dt(1))) - 0.5 .* V .* sqrt(dt(1));


ZC_call_price = discount_alfa* (discount_Talfa_Ti.* normcdf(d1) - strike_i .* normcdf(d2));

end