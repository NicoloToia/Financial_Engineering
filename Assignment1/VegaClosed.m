function vega = VegaClosed(F0,K,KI,B,T,sigma)
% compute the vega of a European Knock-In Call option
%
% INPUT
% F0:       forward price
% K:        strike price
% KI:
% B:        barrier
% T:        time to maturity
% sigma:    volatility
%

% set the change in volatility (1%)
dSigma = 0.01;

if (KI < K) 
    % the option is equivalent to a vanilla call of strike K
    d1 = (log(F0/K) + 0.5*sigma^2*T) / (sigma*sqrt(T));
    vega = B * F0 * sqrt(T) * exp(-0.5*d1^2) / sqrt(2*pi) * dSigma;
else 
    % compute the vega based on the composition of a vanilla with strike KI
    % and a money or nothing with strike KI and value (KI - K)
    d1 = (log(F0/KI) + 0.5*sigma^2*T) / (sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    % Vanilla call with strike KI
    vegaVanillaCall = B * F0 * sqrt(T) * exp(-0.5*d1^2) / sqrt(2*pi) * dSigma;
    % Money or Nothing with strike K and value (KI - K)
    vegaMoneyOrNothing = B * (KI - K) * exp(-0.5*d2^2) / sqrt(2*pi) * dSigma;
    % sum of the two vegas
    vega = vegaVanillaCall + vegaMoneyOrNothing;
end

end