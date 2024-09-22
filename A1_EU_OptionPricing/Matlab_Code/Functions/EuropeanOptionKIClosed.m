 function optionPrice = EuropeanOptionKIClosed(F0,K,KI,B,T,sigma)
% EuropeanOptionKIClosed computes the price of a European option with a
% knock-in barrier using the closed-form solution of Black model
%
% INPUT
% F0:    forward price
% K:     strike
% KI:    barrier
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
%
% OUTPUT
% optionPrice : Price of the option with a knock-in barrier

% if Barrier is less than strike
if (KI < K)
    % the option is equivalent to a Vanilla call option strike K
    optionPrice = EuropeanOptionClosed(F0,K,B,T,sigma, 1);
else % if KI > K
    % Compute the option price by deriving the closed formula by computing
    % under Black's dynamics the price of the composition of a Vanilla call
    % option strike KI and a digital option (cash-or-nothing) with pay off
    % (KI - K)*I(Ftt >= KI)
    d1 = ( log(F0/KI)/(sigma*sqrt(T)) ) + 0.5*sigma*sqrt(T);
    d2 = d1 - sigma*sqrt(T);
    optionPrice = B*(F0*normcdf(d1) - K*normcdf(d2));
end

end