function optionPrice=EuropeanOptionMC(F0,K,B,T,sigma,N,flag)
%European option price with MonteCarlo Simulations
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number of MC simulations
% flag:  1 call, -1 put
%
% OUTPUT
% optionPrice : Price of the option

% Monte Carlo simulation
% Compute the value of the forward at time T for each simulation
Ftt = simulationMC(F0,T,sigma,N);

% Compute the pay-off function for each simulation
payoff = max((Ftt - K),0);

% Compute the price of the call option as the mean of the simulations discounted
% Recall the explicit formula: CallPrice = (1/N)*sum(phi) * B
CallPrice = mean(payoff) * B;

% Exploit put-call parity: C0 - P0 = B*(F0 - K)
if flag == 1
    optionPrice = CallPrice; % Call Price
else 
    optionPrice = CallPrice - B*(F0 - K); % Put Price
end

end % function EuropeanOptionMC