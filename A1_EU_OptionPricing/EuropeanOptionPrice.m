function OptionPrice=EuropeanOptionPrice(F0,K,B,T,sigma,pricingMode,N,flag)
% Option Price with different pricing methods
%
% INPUT:
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% pricingMode: 1 ClosedFormula, 2 CRR, 3 Monte Carlo
% N:     either number of time steps (knots for CRR tree)
%        or number of simulations in MC   
% flag:  1 call, -1 put

if (nargin < 7)
 N = 10000; % Default: N
end 

if (nargin < 8)
 flag = 1; % Default: Call price
end 

switch (pricingMode)
    case 1  % Closed Formula
        OptionPrice = EuropeanOptionClosed(F0,K,B,T,sigma,flag);
    case 2  % CRR
        OptionPrice = EuropeanOptionCRR(F0,K,B,T,sigma,N,flag);
    case 3  % Monte Carlo
        OptionPrice = EuropeanOptionMC(F0,K,B,T,sigma,N,flag);
    otherwise
end
return
end % function EuropeanOptionPrice