function ZC_call_price = ZC_call(strike, fwd_DF, DF, exercise_date, t0, sigma, a, payment_date)

ACT_360 = 2;
ACT_365 = 3;

% compute the volatility using the integral
tau = yearfrac(payment_date, exercise_date, ACT_365);
ttm = yearfrac(t0, exercise_date, ACT_365);

% compute the volatility integral
v_2 = @(t) (sigmaHJM(a, sigma, ttm+tau, t) - sigmaHJM(a, sigma, ttm, t)).^2;

% compute the integral
V = sqrt(1 / ttm  * integral(v_2, 0, ttm));

% define d1 and d2 for every T_i given
d1 = log(fwd_DF / strike) / (V * sqrt(ttm)) + 0.5 * V * sqrt(ttm);
d2 = d1 - V * sqrt(ttm);

ZC_call_price = DF * (fwd_DF * normcdf(d1) - strike * normcdf(d2));

end