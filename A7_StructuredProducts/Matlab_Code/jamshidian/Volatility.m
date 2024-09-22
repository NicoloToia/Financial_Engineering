function vol = Volatility(t_alpha,t0, Ti, sigma, a)

ACT_365 = 3;

Ti = yearfrac(t0,Ti,ACT_365);
talfa = yearfrac(t0,Talfa,ACT_365);

integrand = @(u) (sigma/a * (1 - exp(-a * (Ti-u)))-sigma/a * (1 - exp(-a * (talfa-u))) ).^2;

int = integral(integrand, 0, talfa,"ArrayValued",true);

vol = sqrt(int./talfa);

end