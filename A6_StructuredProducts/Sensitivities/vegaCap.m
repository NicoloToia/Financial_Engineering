function vega = vegaCap(strike, cap_ttm, spot_vols, spot_ttms, mkt_vols, ttms, strikes, discounts, dates)
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

% compute the dates for the cap
exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
    calmonths(3:3:cap_ttm*12-3)';
payment_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
    calmonths(6:3:cap_ttm*12)';

% move to business days if needed
exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
    busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
% move to business days if needed
payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
    busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());

% conver to datenum
exercise_dates = datenum(exercise_dates);
payment_dates = datenum(payment_dates);

% compute the cap price
cap_price_0 = CapSpot(strike, exercise_dates, payment_dates, spot_vols, spot_ttms, strikes, ...
    discounts, dates);

% shift the spot vols by 1 bp in the corresponding row
shift = 0.0001;

% find the row corresponding to the cap_ttm
row = find(ttms == cap_ttm);

% shift the row by 1 bp
shifted_vols = mkt_vols;
shifted_vols(row, :) = shifted_vols(row, :) + shift;

% recompute the mkts prices
mkt_prices = MarketCapPrices(ttms, strikes, shifted_vols, discounts, dates);

% recalibrate the spot vols
[shifted_ttms, shifted_vols] = spotVols(mkt_prices, ttms, strikes, shifted_vols, discounts, dates);

% compute the cap price with the shifted vols
cap_price_shift = CapSpot(strike, exercise_dates, payment_dates, shifted_vols, shifted_ttms, strikes, ...
    discounts, dates);

% compute the vega
vega = (cap_price_shift - cap_price_0);

end