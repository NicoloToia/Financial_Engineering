function optionPrice=EuropeanOptionClosed(F0,K,B,T,sigma,flag)
%European option price with Closed formula by Black model
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% flag:  1 call, -1 put
%
% OUTPUT
% optionPrice : Price of the option

[call, put] = blkprice(F0, K, 0, T, sigma);

if flag == 1 
   optionPrice= B*call;
else
   optionPrice= B*put;
end

end % function EuropeanOptionClosed