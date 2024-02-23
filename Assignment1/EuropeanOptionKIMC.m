function optionPrice=EuropeanOptionKIMC(F0,K,KI,B,T,sigma,N)
% EuropeanOptionKIMC computes the price of a European option with a
% knock-in barrier using Monte Carlo method
%
% INPUT
% F0 :   Initial price of the underlying asset
% K :    Strike price
% KI :   Barrier level
% B :    Discount factor
% T :    Time to expiration
% sigma: Volatility
% N :    Number of MC simulations
%
% OUTPUT
% optionPrice : Price of the option with a knock-in barrier

% Compute the value of the forward by Monte Carlo simulation
Ftt = simulationMC(F0,T,sigma,N);

% Compute the prices by discounting the mean of the payoff
payoff = max((Ftt - K), 0) .* (Ftt > KI);
optionPrice = B * mean(payoff);

end