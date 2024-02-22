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
    d1 = (log(F0/KI) + 0.5*sigma^2*T)/(sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    optionPrice = B*(F0*normcdf(d1) - K*normcdf(d2));
end

end