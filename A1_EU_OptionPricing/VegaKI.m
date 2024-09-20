function vega=VegaKI(F0,K,KI,B,T,sigma,N,flagNum)
% Compute the vega of a knock-in call option
%
% INPUT:
% F0:      forward price
% K:       strike price
% KI:      knock-in barrier
% B:       discount factor
% T:       maturity
% sigma:   volatility
% N:       number of simulations for MC and number of steps for CRR
% flagNum: 1 for CRR, 2 for MC, 3 for closed-form
%
% OUTPUT:
% vega:    vega of the option

% Redirect to specic functions
switch (flagNum)
    case 1
        vega = VegaCRR(F0,K,KI,B,T,sigma,N);
    case 2
        vega = VegaMC(F0,K,KI,B,T,sigma,N);
    case 3
        vega = VegaClosed(F0,K,KI,B,T,sigma);
    otherwise
        error('VegaKI: invalid flagNum');
end

end