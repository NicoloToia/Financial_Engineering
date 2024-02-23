function vega = VegaClosed(F0,K,KI,B,T,sigma)
% Compute the vega of a European Knock-In Call option by the closed formula
%
% INPUT
% F0:    forward price
% K:     strike price
% KI:    knock-in barrier
% B:     barrier
% T:     time to maturity
% sigma: volatility
%
% OUTPUT:
% vega:  vega of the option

% Set the change in volatility at (1%)
dSigma = 0.01;

if (KI < K) 
    % The option is equivalent to a vanilla call of strike K
    d1 = (log(F0/K) + 0.5*sigma^2*T) / (sigma*sqrt(T));
    vega = B * F0 * sqrt(T) * exp(-0.5*d1^2) / sqrt(2*pi) * dSigma;
else % if KI > K
    % Compute the vega based on the composition of a vanilla with strike KI
    % and a money or nothing with strike KI and value (KI - K)
    d1 = (log(F0/KI) + 0.5*sigma^2*T) / (sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);

%    vega = B * F0 * sqrt(T) * exp(-d1^2/2) / sqrt(2*pi);
    der_d1 = -1/sigma^2 * log(F0/KI)/sqrt(T) + 0.5*sqrt(T); % derivative of d1
    der_d2 = -1/sigma^2 * log(F0/KI)/sqrt(T) - 0.5*sqrt(T); % derivative of d2

    vega = B * ( F0 * (exp(-d1^2/2) / sqrt(2*pi) ) * der_d1 - K * (exp(-d2^2/2) / sqrt(2*pi) ) * der_d2 ) * dSigma;
end

end