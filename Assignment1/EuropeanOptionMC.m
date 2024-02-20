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

% generate random numbers from a standard Gaussian distribution
g = randn(N,1);

priceMC = 0;
for i = 1:N
    priceMC = priceMC + (B/N)*max(F0*exp(-(sigma^2)*T*0.5 + sigma*sqrt(T)*g(i)) - K,0); 
end 
if flag == 1
    optionPrice = priceMC;
else
    % Call-Put parity (European)
    optionPrice = priceMC - B*(F0 - K);
end

end % function EuropeanOptionMC
