function [M,stdEstim]=PlotErrorMCAV(F0,K,B,T,sigma)
% Error plot for MCAV method
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
%
% OUTPUT
% M :       number of simulations (1x20 vector)
% stdEstim: MC price standard deviation estimate for each M (1x20)

% Compute the number of simulations to use
m=1:20;
M=2.^m;

% Inizialize stdEstim vector
stdEstim=zeros(1,20);

for i=1:length(M)

    % Monte Carlo simulation AV
    [Ftt, FttAV] = simulationMCAV(F0,T,sigma,M(i));
    
    % Compute the payoff
    payoff =  0.5 * (max(Ftt-K,0) + max(FttAV-K,0));

    % Compute the variance of the price
    stdEstim(i) = B * std(payoff) / sqrt(M(i));
end

end