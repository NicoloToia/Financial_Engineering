function optionPrice = EuropeanOptionMCAV(F0,K,B,T,sigma,N,flag)
% EuropeanOptionMCAV: European option price using Monte Carlo simulation
% with antithetic variates
%
% INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations
% flag:  1 call, -1 put
%
% OUTPUT
% optionPrice : Price of the option

% Monte Carlo simulation with antithetic variates
[Ftt, FttAV] = simulationMCAV(F0,T,sigma,N);

% compute the discounted payoffs
discPayoff = B*max(flag*(Ftt-K),0);
discPayoffAV = B*max(flag*(FttAV-K),0);

% sum and average the discounted payoffs
optionPrice = mean(0.5 * (discPayoff + discPayoffAV));

end