function optionPrice = EuropeanOptionKIClosed(F0,K,KI,B,T,sigma)

% EuropeanOptionKIClosed computes the price of a European option with a
% knock-in barrier using the closed-form solution of Black model
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% KI:    barrier
% T:     time-to-maturity
% sigma: volatility
% N:     either number of time steps (knots for CRR tree)

% if Barrier is less than strike
if (KI < K)
    % the option is equivalent to a Vanilla call option
    optionPrice = EuropeanOptionClosed(F0,K,B,T,sigma, 1);
else % if KI > K
    % Option is equivalent to a Vanilla call option of strike KI
    % plus a Cash-or-nothing call option of strike KI and cash amount KI-K
    d2 = @(F0,K,T,sigma) (log(F0/K) - 0.5*sigma^2*T)/(sigma*sqrt(T));
    optionPrice = EuropeanOptionClosed(F0,KI,B,T,sigma, 1) + (KI-K) * normcdf(d2(F0,KI,T,sigma));
end

% discounting
optionPrice = B * optionPrice;

end