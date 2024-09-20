function swapRateHW(r, fwd_DF, sigma, a, t, t0)

sigma_hjm = sigma/a * (1 - exp(-a * dt)) % to be modified in dt

integral_ZC = (sigma/a)^2 (t - t0 - 3/(2*a) + 1/a * exp(-a(t-t0)) * (2 - 0.5*exp(-a*(t -t0))) )

% compute the ZCBs
B = fwd_DF * exp( - r * )

end