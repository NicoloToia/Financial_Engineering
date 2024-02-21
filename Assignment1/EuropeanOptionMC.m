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

% % Monte Carlo simulation (one time step)
% % Compute the value of the forward at time T for each simulation
Ftt=simulationMC(F0,T,sigma,N);

% Compute the call option price for each simulation
Ct = max((Ftt - K),0);

% Compute the price of the call option as the mean of the simulations discounted
C0 = mean(Ct) * B;


if flag == 1
    optionPrice = C0;
else % put option
    % leverage the put-call parity
    % C0 - P0 = B*(F0 - K)
    optionPrice = C0 - B*(F0 - K);
end

end % function EuropeanOptionMC
