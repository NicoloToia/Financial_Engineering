function vega = VegaMC(F0,K,KI,B,T,sigma,N);
% compute the vega of KI European call option using Monte Carlo simulation
%
% INPUT
% F0: initial stock price
% K: strike price
% KI: knock-in barrier
% B; discount factor
% T: time to maturity
% sigma: volatility
% N: number of simulations

% use center scheme to compute the vega
dSigma = 0.01;
prevPrice = EuropeanOptionKIMC(F0,K,KI,B,T,sigma-dSigma,N);
nextPrice = EuropeanOptionKIMC(F0,K,KI,B,T,sigma+dSigma,N);

% f'(x) = (f(x+h) - f(x-h)) / (2h)
vega = (nextPrice - prevPrice) / (2*dSigma);

end