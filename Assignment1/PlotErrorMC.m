function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% error plot for CRR method
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility

m=[1:20];
M=2.^m;
stdEstim=zeros(1,20);

for i=1:length(M)
    g = randn(M(i),1);
    % Monte Carlo simulation (one time step)
    % Compute the value of the forward at time T for each simulation
    % Black Model: Ft = F0 * exp(-(sigma^2)*T*0.5 + sigma*sqrt(T)*g)
    Ftt = F0 * exp( -0.5 * sigma^2 * T  + sigma * sqrt(T) * g);
    
    % Compute the call option price for each simulation
    stdEstim(i)=B^2*std(max((Ftt - K),0))/sqrt(M(i));
end
figure
title('MC')
loglog(M,stdEstim)
hold on
loglog(M,1./sqrt(M))
end