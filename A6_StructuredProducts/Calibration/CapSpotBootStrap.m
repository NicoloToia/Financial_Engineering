function [Price, sigmas] = CapSpotBootStrap(Strike, sigma_alpha, T_alpha, sigma_beta, caplet_ttms, caplet_yf, caplet_DF, Libor)
% CapFlat: Compute the price of a cap using the flat volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%)
%   sigma_alpha : Volatility of the previous caplet
%   sigma_beta : Volatility of final caplet
%   T_alpha : Date of the previous caplet
%   caplet_ttms : Caplet maturities
%   caplet_yf : Caplet year fractions
%   caplet_DF : Caplet discount factors
%   Libor : forward Libor rates

% linearize the volatilities
sigmas = sigma_alpha + (sigma_beta - sigma_alpha) / (caplet_ttms(end) - T_alpha) ...
    * (caplet_ttms - T_alpha);

% apply the bachelier formula in the vector form
K = Strike / 100 * ones(size(caplet_ttms));
dn = (Libor - K) ./ (sigmas .* sqrt(caplet_ttms));

term_1 = (Libor - K) .* normcdf(dn);
term_2 = sigmas .* sqrt(caplet_ttms) .* normpdf(dn);

% compute the price of the cap
Price = sum(caplet_DF .* caplet_yf .* (term_1 + term_2));

end