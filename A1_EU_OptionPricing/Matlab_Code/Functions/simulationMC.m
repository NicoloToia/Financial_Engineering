function Ftt=simulationMC(F0,T,sigma,N)
% MonteCarlo simulations forward
%
% INPUT 
% F0:    forward price
% T:     time-to-maturity
% sigma: volatility
% N:     number of MC simulations
%
% OUTPUT
% Ftt :  simulated forward at time T

% Extract N random numbers from a standard gaussian
g = randn(N,1);

% Monte Carlo simulation 
% Compute the value of the forward at time T
% Black Model: Ft = F0 * exp(-(sigma^2)*T*0.5 + sigma*sqrt(T)*g)
Ftt = F0 * exp( -0.5 * sigma^2 * T  + sigma * sqrt(T) * g);
end