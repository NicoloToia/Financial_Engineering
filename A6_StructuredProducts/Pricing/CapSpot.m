function Price = CapSpot(strike, caplet_ttms, caplet_yf, caplet_DF, Libor, spot_vols, ttms, strikes)
% CapSpot: Compute the price of a cap using the spot volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%) 
%   caplet_ttms : Caplet maturities
%   caplet_yf : Caplet year fractions
%   caplet_DF : Caplet discount factors
%   Libor : forward Libor rates
%   spot_vols : Spot volatilities (calibrated to the market caplets)
%   ttms : Time to maturities of the caplets (in years)
%   strikes : Strikes of the caplets


% interpolate the spot volatilities to the target strikes and deltas
sigmas = intSpotVols(strike, caplet_ttms, spot_vols, ttms, strikes);

% apply the bachelier formula in the vector form
K = strike / 100 * ones(size(caplet_ttms));

dn = (Libor - K) ./ (sigmas .* sqrt(caplet_ttms));

term_1 = (Libor - K) .* normcdf(dn);
term_2 = sigmas .* sqrt(caplet_ttms) .* normpdf(dn);

% compute the price of the cap
Price = sum(caplet_DF .* caplet_yf .* (term_1 + term_2));

end