function optionPrice = EuropeanOptionKIMC(F0,K,KI,B,T,sigma,N)

% EuropeanOptionKIMC: Up & In European call option
%
% INPUT
% F0 : Initial price of the underlying asset
% K : Strike price
% KI : Barrier level
% B : Discount factor
% T : Time to expiration
% sigma : Volatility
% N : Number of simulations
%
% OUTPUT
% optionPrice : Price of the option

    % simulate the forward price
    Ftt=simulationMC(F0,T,sigma,N);

    % compute the prices
    payoff = max((Ftt - K), 0) .* (Ftt > KI);
    optionPrice = B * mean(payoff);

end