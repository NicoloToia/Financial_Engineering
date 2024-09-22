function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)
% Error plot for MC method. As error an estimation of the unbiased standard
% deviation of the MC price is analyzed
%
%INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
%
% OUTPUT
% M :       number of simulations (1x20 vector)
% stdEstim: MC price standard deviation estimate for each M (1x20)

% Compute the number of steps to use
m=1:20;
M=2.^m;

% Inizialize stdEstim vector
stdEstim=zeros(1,20);

for i=1:length(M)
    % Monte Carlo simulation (one time step) for each M(i)
    Ftt = simulationMC(F0, T, sigma, M(i));
    
    % Compute the variance of the price
    stdEstim(i) = B * std( max((Ftt - K),0) ) / sqrt(M(i));
end

% spread is 1 bp
spread = 10^-4;

% Plot the results of MC
figure
loglog(M,stdEstim)
title('MC Error')
xlabel('M'); ylabel('errorMC')
hold on 
loglog(M,1./sqrt(M))
% cutoff is based on the spread
loglog(M, spread * ones(length(M),1))
legend('MC','1/sqrt(M)','cutoff')

end