function vega = VegaMCestim(F0,K,KI,B,T,sigma,N)
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

dSigma = 0.01;

% simulation the forward
Ftt = simulationMC(F0,T,sigma,N);

% compute the estimator
estimator = (B + KI - K) * ((log(Ftt./F0) - 0.5*sigma^2*T) / sigma) .* Ftt .* (Ftt>KI);

% compute the vega
vega = mean(estimator) * dSigma;

end