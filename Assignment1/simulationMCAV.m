function [Ftt, FttAV]=simulationMCAV(F0,T,sigma,N)
% MonteCarlo simulations forward with Antithetic Variables
%
%INPUT 
% F0:    forward price
% T:     time-to-maturity
% sigma: volatility
% N:     number of simulations

% extract N random numbers from a standard gaussian
g = randn(N,1);

Ftt = F0*exp(-0.5*sigma^2*T+sigma*sqrt(T)*g);
FttAV = F0*exp(-0.5*sigma^2*T-sigma*sqrt(T)*g);

end