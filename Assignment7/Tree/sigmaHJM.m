function sigmaHJM = sigmaHJM(a, sigma, t, s)
% compute the HJM volatility

sigmaHJM = sigma/a * (1 - exp(-a * (t - s)));

end