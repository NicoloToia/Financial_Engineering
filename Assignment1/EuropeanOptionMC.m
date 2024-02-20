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
%S0 = F0/(exp(-d*TTM)/B);
S0 = 1;
if flag == 1

    S(:) = S0*exp(g*sigma*sqrt(T)); 
else
    % put call parity
end

optionPrice = B*mean(max(S(1,:) - K,0));

end
