function Vega = total_Vega(mkt_vols, ttms, strikes, spol_A, fixed_rate_B, spol_B, cap_5y, cap_10y, cap_15y, discounts, dates, mkt_prices)

    % shock the volatility of the caps by 1bp
    bp = 0.0001;
    shock = bp * ones(size(mkt_vols));
    shockVols = mkt_vols + shock;
    spot_vols = spotVols(mkt_prices, ttms, strikes, mkt_vols, discounts, dates);
    price_0 = computeUpfront(spot_vols, ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    price_shift = computeUpfront(shockVols, ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
            cap_5y, cap_10y, cap_15y, discounts, dates);
    % compute the total vega
    Vega = (price_shift - price_0)/price_0;
end