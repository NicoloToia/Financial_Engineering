function sigmaHJM = sigmaHJM(a, sigma, t, s)
% compute the HJM volatility

ACT_365 = 3;

sigmaHJM = sigma/a * (1 - exp(-a * yearfrac(s,t, ACT_365)));

end