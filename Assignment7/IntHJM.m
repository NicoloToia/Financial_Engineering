function IntHJM = IntHJM(a, sigma, t_i, t0)

IntHJM = (sigma/a)^2 * (t_i - t0 - 3/(2*a) + 1/a * exp(-a*(t_i-t0)) * (2 - 0.5*exp(-a*(t_i -t0))) );

end