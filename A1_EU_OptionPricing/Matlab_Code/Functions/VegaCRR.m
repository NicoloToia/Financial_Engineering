function vega = VegaCRR(F0,K,KI,B,T,sigma,N)
% Compute the vega of KI European call option using CRR Method
%
% INPUT
% F0:    initial stock price
% K:     strike price
% KI:    knock-in barrier
% B:     discount factor
% T:     time to maturity
% sigma: volatility
% N:     number of time steps
%
% OUTPUT:
% vega:  vega of the option

% Set the change in volatility at (1%)
dSigma = 0.01;

% Use center scheme to compute the vega
prevPrice = EuropeanOptionKICRR(F0,K,KI,B,T,sigma-dSigma,N);
nextPrice = EuropeanOptionKICRR(F0,K,KI,B,T,sigma+dSigma,N);

% f'(x) = (f(x+h) - f(x-h)) / (2h)
vega = ((nextPrice - prevPrice) / (2*dSigma)) * dSigma;

end