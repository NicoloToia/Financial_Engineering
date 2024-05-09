function IntHJM = IntHJM(a, sigma, t_i, t0, tau)

% define the integrand function
integrand = @(u) ( (sigmaHJM(a, sigma, (t_i + tau), u)).^2 - (sigmaHJM(a, sigma, t_i, u)).^2 );

% compute the integral
IntHJM = integral(integrand, t0, t_i);

% IntHJM = (sigma/a)^2 * (t_i - t0 - 3/(2*a) + 1/a * exp(-a*(t_i-t0)) * (2 - 0.5*exp(-a*(t_i -t0))) );

end