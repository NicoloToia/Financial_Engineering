function sensitivities = vegaBuckets(mkt_vols, ttms, strikes, X_0, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, caplet_ttms, caplet_yf, caplet_DF, Libor)
% VEGABUCKETS computes the vega bucket sensitivities for the certificate
%
% INPUTS
%   mkt_vols: market volatilities
%   ttms: time to maturities of the mkt_vols
%   strikes: strikes
%   X_0: upfront payment
%   spol_A: spread over libor for the first quarter
%   fixed_rate_B: fixed rate for the second quarter
%   spol_B: spread over libor for the third quarter
%   cap_5y: strike of the cap from 0 to 5y
%   cap_10y: strike of the cap from 5 to 10y
%   cap_15y: strike of the cap from 10 to 15y
%   discounts: discounts
%   dates: dates of the market data

% initialize the sensitivities
sensitivities = zeros(length(mkt_vols), 1);
% shift is 1 bp
shift = 10^(-4);
% select only the needed data for the certificate
cf_caplet_ttms = caplet_ttms(1:15*4);
cf_caplet_yf = caplet_yf(1:15*4);
cf_caplet_DF = caplet_DF(1:15*4);
cf_libor = Libor(1:15*4);

% for each maturity, compute the vega bucket sensitivity
for i = 1:length(mkt_vols)
    
    % shift the vols of 1 bp for the corresponding row
    shifted_vols = mkt_vols;
    shifted_vols(i, :) = shifted_vols(i, :) + shift;

    % recompute the market prices
    mkt_prices = MarketCapPrices(ttms, strikes, shifted_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);

    % recalibrate the spot volatilities
    shift_vols = spotVols(mkt_prices, ttms, strikes, shifted_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);

    % recompute the upfront payment
    X_shift = computeUpfront(shift_vols, caplet_ttms, strikes, spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, cf_caplet_ttms, cf_caplet_yf, cf_caplet_DF, cf_libor);
    % compute the vega bucket sensitivity
    sensitivities(i) = (X_shift - X_0);
end

end