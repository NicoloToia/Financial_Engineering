function [Call, Put] = blackPrice(F_0, K, T, sigma, B_0, flag)
% Compute the price of a call and put option using the Black-Scholes formula
%
% INPUT:
%   F_0: forward price at time 0
%   K: strike price
%   T: time to maturity
%   sigma: volatility
%   B_0: discount factor at time 0
%   flag: flag to choose the option type ("call" or "put")
%
% OUTPUT:
%   Call: price of the call option
%   Put: price of the put option

% compute the d_1 and d_2
d_1 = (log(F_0 / K) + (0.5 * sigma^2) * T) / (sigma * sqrt(T));
d_2 = d_1 - sigma * sqrt(T);

% compute the call and put prices
Call = B_0 * (F_0 * normcdf(d_1) - K * normcdf(d_2));
Put = B_0 * (K * normcdf(-d_2) - F_0 * normcdf(-d_1));

end