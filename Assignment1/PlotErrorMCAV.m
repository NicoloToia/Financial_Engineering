function [M,stdEstim]=PlotErrorMCAV(F0,K,B,T,sigma)
% error plot for CRR method
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

    % monte carlo simulation AV
    [Ftt, FttAV] = simulationMCAV(F0,T,sigma,M(i));
    
    % compute the payoff
    payoff =  0.5 * (max(Ftt-K,0) + max(FttAV-K,0));

    stdEstim(i) = B * std(payoff) / sqrt(M(i));
end

end