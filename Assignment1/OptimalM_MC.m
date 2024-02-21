function [M_optimal]=OptimalM_MC(F0,K,B,T,sigma, bidask)
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
stdEstim=zeros(1,length(m));
for i=1:length(M)
    g = randn(M(i),1);
    % Monte Carlo simulation (one time step)
    % Compute the value of the forward at time T for each simulation
    % Black Model: Ft = F0 * exp(-(sigma^2)*T*0.5 + sigma*sqrt(T)*g)
    Ftt = F0 * exp( -0.5 * sigma^2 * T  + sigma * sqrt(T) * g);
    
    % Compute the call option price for each simulation
    stdEstim(i)=B*std(max((Ftt - K),0))/sqrt(M(i)); 
end
ii = 1; 
while (stdEstim(ii)>bidask && ii<21)
    ii = ii+1; 
end 
M_optimal = ii; 
end