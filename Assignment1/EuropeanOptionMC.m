function optionPrice=EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
%European option price with MonteCarlo Simulations
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations
% flag:  1 call, -1 put

% generate random numbers
g = randn(N,1);
% compute St based on the forward dynamics
St = F0 * exp( -0.5 * sigma^2 * T + sigma * sqrt(T) * g);

% payoff of the option
Payoff = max(flag*(St-K), 0);

% compute the option price
optionPrice = B * mean(Payoff);

end
