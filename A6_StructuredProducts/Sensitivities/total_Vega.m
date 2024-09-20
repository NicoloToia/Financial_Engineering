function Vega = total_Vega(mkt_vols, ttms, strikes, X_0, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, caplet_ttms, caplet_yf, caplet_DF, Libor)
% TOTAL_VEGA computes the total vega of the certificate
%
% INPUTS
%   mkt_vols: market volatilities
%   ttms: time to maturities
%   strikes: strikes
%   spol_A: spread over libor for the first quarter
%   fixed_rate_B: fixed rate for the second quarter
%   spol_B: spread over libor for the third quarter
%   cap_5y: strike of the cap from 0 to 5y
%   cap_10y: strike of the cap from 5 to 10y
%   cap_15y: strike of the cap from 10 to 15y
%   caplet_ttms: caplet maturities
%   caplet_yf: caplet year fractions
%   caplet_DF: caplet discount factors
%   Libor: forward Libor rates

% shock the volatility of the caps by 1bp
shift = 10^(-4);
% shift the market vols by 1bp
shift_vols = mkt_vols + shift;

% recompute the market prices
mkt_prices = MarketCapPrices(ttms, strikes, shift_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);

% recalibrate the spot volatilities
spot_vols = spotVols(mkt_prices, ttms, strikes, shift_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);
% select only the needed data for the certificate
cf_caplet_ttms = caplet_ttms(1:15*4);
cf_caplet_yf = caplet_yf(1:15*4);
cf_caplet_DF = caplet_DF(1:15*4);
cf_libor = Libor(1:15*4);
% recalculate the upfront
X_1 = computeUpfront(spot_vols, caplet_ttms, strikes, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, cf_caplet_ttms, cf_caplet_yf, cf_caplet_DF, cf_libor);
% compute the total vega
Vega = (X_1 - X_0);

end