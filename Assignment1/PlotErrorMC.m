function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% error plot for MC method
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility

m=1:20;
M=2.^m;
stdEstim=zeros(1,20);

for i=1:length(M)
    % Monte Carlo simulation (one time step)
    Ftt = simulationMC(F0, T, sigma, M(i));
    
    % Compute the variance of the price
    stdEstim(i) = B * std( max((Ftt - K),0) ) / sqrt(M(i));
end

end