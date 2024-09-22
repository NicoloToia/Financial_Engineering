function sensitivities = vegaBucketsMatrix(mkt_vols, ttms, strikes, X_0, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, discounts, dates)
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
sensitivities = zeros(size(mkt_vols));
% shift is 1 bp
shift = 0.0001;

% for each volatility, compute the vega bucket sensitivity
for i = 1:length(mkt_vols)
    for j = 1:width(mkt_vols)
        % shift the volatility of 1 bp
        shifted_vols = mkt_vols;
        shifted_vols(i, j) = shifted_vols(i, j) + shift;
        % recompute the cap prices
        mkt_prices = MarketCapPrices(ttms, strikes, shifted_vols, discounts, dates);
        % recalibrate the spot vols
        [shifted_ttms, shifted_vols] = spotVols(mkt_prices, ttms, strikes, shifted_vols, discounts, dates);
        % recompute the upfront payment
        X_shift = computeUpfront(shifted_vols, shifted_ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
            cap_5y, cap_10y, cap_15y, discounts, dates);
        % compute the vega bucket sensitivity
        sensitivities(i, j) = (X_shift - X_0);
    end
end

end