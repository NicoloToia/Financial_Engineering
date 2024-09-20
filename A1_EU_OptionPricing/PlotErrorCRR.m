function [M, errorCRR]=PlotErrorCRR(F0,K,B,T,sigma)
% Error plot for CRR method.  As error is considered the difference 
% in absolute value w.r.t. to the exact value.
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
%
% OUTPUT
% M :       number of steps (1x10 vector)
% errorCRR: error for each number of steps (1x10)

% Compute the price of a call option using Black's formula and fix it as
% exact price
% M does not matter, the closed formula is invariant
exactPrice = EuropeanOptionPrice(F0,K,B,T,sigma,1,0,1);

% Compute the number of steps to use
m = 1:10;
M = 2.^m;

% Compute the error for each number of steps (row vector)
errorCRR = zeros(1, length(M));

for i = 1:length(M)
    % compute the price of a call using the CRR method for M(i)
    priceCRR = EuropeanOptionPrice(F0,K,B,T,sigma,2,M(i),1);
    % compute the error as the absolute vlaue of the difference
    errorCRR(i) = abs(priceCRR - exactPrice);
end

% spread is 1 bp
spread = 10^-4;

% Plot the results of CRR
figure
loglog(M,errorCRR)
title('CRR Error')
xlabel('M'); ylabel('errorCRR')
hold on
loglog(M, 1./M)
% cutoff is based on the spread
loglog(M, spread * ones(length(M),1))
legend('CRR','1/M','cutoff')

end