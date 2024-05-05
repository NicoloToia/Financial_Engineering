function vega = totalVegaCap(cap_0_5y, cap_ttms, cap_yf, cap_DF, cap_libor, mkt_vols, ttms, strikes, ...
    caplet_ttms, caplet_yf, caplet_DF, Libor)
% VEGACAP computes the vega of a cap from 0 to cap_ttm with given strike
%
% INPUTS
%   strike_5y: strike of the cap from 0 to 5y
%   cap_ttm: time to maturity of the cap (in years)
%   discounts: discounts
%   dates: dates of the market data
%   spot_vols: spot volatilities (calibrated from the market data)
%   spot_ttms: time to maturities of the spot_vols
%   strikes: strikes
%   mkt_vols: market volatilities (to be shifted)

% shift the spot vols by 1 bp in the corresponding row
shift = 10^(-4);

% shift the row by 1 bp
shifted_mkt_vols = mkt_vols + shift;

% recompute the mkts prices
mkt_prices = MarketCapPrices(ttms, strikes, shifted_mkt_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);

% recalibrate the spot vols
shifted_spot_vols = spotVols(mkt_prices, ttms, strikes, shifted_mkt_vols, caplet_ttms, caplet_yf, caplet_DF, Libor);

% compute the cap price with the shifted vols
cap_price_shift = CapSpot(strike, cap_ttms, cap_yf, cap_DF, cap_libor, shifted_spot_vols, caplet_ttms, strikes);

% compute the vega
vega = (cap_price_shift - cap_price_0);

end