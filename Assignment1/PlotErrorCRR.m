function [M, errorCRR] = PlotErrorCRR(F0,K,B,T,sigma)
% error plot for CRR method
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
%
%OUTPUT
% M:    number of steps (1x10 vector)
% errorCRR: error for each number of steps (1x10)

% compute the price of a call option using Black's formula
% M does not matter, the closed formula is invariant
exactPrice = EuropeanOptionPrice(F0,K,B,T,sigma,1,0,1);

% compute the number of steps to use
m = 1:10;
M = 2.^m;

% compute the error for each number of steps (row vector)
errorCRR = zeros(1, length(M));

for i = 1:length(M)
    % compute the price of a call using the CRR method for M
    priceCRR = EuropeanOptionPrice(F0,K,B,T,sigma,2,M(i),1);
    % compute the error as the absolute vlaue
    errorCRR(i) = abs(priceCRR - exactPrice);
end

% spread is 1 bp
spread = 10^-4;

% Plot the results of CRR
figure
subplot(1,2,1)
loglog(M,errorCRR)
title('CRR')
hold on
loglog(M, 1./M)
% cutoff
loglog(M, spread * ones(length(M),1))
legend('CRR','1/nCRR','cutoff')

end