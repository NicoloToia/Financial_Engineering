% runAssignment6_Group16
%  group 16, AY2023-2024
% 
%
% to run:
% > runAssignment6_Group16

% clear workspace
clear all;
close all;
clc;

% set the clock to find the time of execution
tic;

%% Add the directories to the path

addpath('Data');
addpath('Bootstrap');
addpath('Pricing')
addpath('Calibration');
addpath('Sensitivities');
addpath('Hedging');
addpath('Plotting');

%% Settings

formatData='dd/mm/yyyy'; % Pay attention to your computer settings 

%% Read market data

% This function works on Windows OS. Pay attention on other OS.
[datesSet, ratesSet] = readRatesData('MktData_CurveBootstrap_20-2-24', formatData);

[ttms, strikes, mkt_vols] = readVolData('Caps_vol_20-2-24');

%% Construct the swap dates
 
swapDates = datetime(datesSet.settlement, 'ConvertFrom', 'datenum') + calyears(1:50)';
swapDates(~isbusday(swapDates, eurCalendar())) = busdate(swapDates(~isbusday(swapDates, eurCalendar())), 'modifiedfollow', eurCalendar());
swapDates = datenum(swapDates);

%% Interpolate the mid market rates for swaps from the market data

% find the actually quoted swap dates
datesSet.swaps = [swapDates(1:12); swapDates(15:5:30); swapDates(40:10:50)];

% save the quoted dates (those used in the bootstrap)
quoted_dates = [datesSet.depos(1:4); datesSet.futures(1:7, 2); datesSet.swaps(2:end)];

%% Bootstrap the discount factors from the market data

% Bootstrap the discount factors from the market data
[dates, discounts] = bootstrap(datesSet, ratesSet);

% compute zero rates from the discounts founded in the bootstrap
zeroRates = zeroRates(dates, discounts);

%% Plot Results

% plot_rates_discounts(dates, discounts, zeroRates);

%% Obtain the Cap Prices from the market data via Bachelier formula

% Cap prices
mkt_cap_prices = MarketCapPrices(ttms, strikes, mkt_vols, discounts, dates);

% plot the cap prices
% plotCaps(mkt_cap_prices, ttms, strikes);

%% Compute the spot volatilites

[spot_ttms, spot_vols] = spotVols(mkt_cap_prices, ttms, strikes, mkt_vols, discounts, dates);

final_vols = spot_vols(4*ttms-1,:);
final_ttms = spot_ttms(4*ttms-1,:);

%% Plot the spot volatilities surface against the flat volatilities

% plot_vols(spot_vols, spot_ttms, mkt_vols, ttms, strikes);

%% Price the certificate

spol_A = 2; % 2% spread for party A
fixed_rate_B = 3; % 3% fixed rate for party B (first quarter only)
spol_B = 1.1; % 1.1% spread for party B
cap_rate_5y = 4.3; % 4.3% cap rate for the coupon of party B from 0 to 5 years
cap_rate_10y = 4.6; % 4.6% cap rate for the coupon of party B from 5 to 10 years
cap_rate_15y = 5.1; % 5.1% cap rate for the coupon of party B from 10 to 15 years

X = computeUpfront(final_vols, final_ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
    cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate ---')
disp(['The upfront payment is: ', num2str(X*100), '%']);
Notional = 50 * 10^6;
disp(['The upfront payment is: ', num2str(X*Notional), ' EUR']);
disp('--- --- ---')

%% Delta-bucket sensitivity

% compute the delta-bucket sensitivity
[delta_dates, delta_buckets] = deltaBuckets(datesSet, ratesSet, quoted_dates, spot_vols, spot_ttms, strikes, X, dates(1), spol_A, ...
    fixed_rate_B, spol_B, cap_rate_5y, cap_rate_10y, cap_rate_15y);

disp('--- Delta-bucket sensitivity of the Certificate ---')
disp('Date | DV01 (Notional = 100) | DV01 (EUR) |')
for i = 1:length(delta_dates)
    disp([datestr(delta_dates(i), formatData), ' | ', num2str(delta_buckets(i) * 100), ' | ', ...
        num2str(delta_buckets(i) * Notional)]);
end
disp('--- --- ---')

%% Plot the delta-bucket sensitivities

% plot_delta_buckets(delta_dates, delta_buckets);
% plot_delta_buckets(delta_dates, delta_buckets * Notional);
% these values are expressed in bp, to obtain EUR values, multiply by the notional and by 1 bp

%% Total Vega

total_vega = total_Vega(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
    cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);
% print the total vega
disp('--- Total vega of the Certificate ---')
disp(['The total vega is: ', num2str(total_vega * 100), ' (Notional = 100) EUR']);
disp(['The total vega is: ', num2str(total_vega * Notional), ' EUR']);
disp('--- --- ---')

%% Vega Buckets sensitivity matrix

if ~isfile('Data/vega_matrix.mat')
    % compute the vega buckets matrix
    vega_matrix = vegaBucketsMatrix(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
        cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);
    save('Data/vega_matrix.mat', 'vega_matrix')
else
    load('Data/vega_matrix.mat');
end

% plot_Vega_matrix(vega_matrix, ttms, strikes)

%% Vega bucket sensitivity

% compute the vega dates from the ttms
vega_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(ttms);
% move to business days if needed
vega_dates(~isbusday(vega_dates, eurCalendar())) = ...
    busdate(vega_dates(~isbusday(vega_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
vega_dates = datenum(vega_dates);

% compute the vega bucket sensitivities
vega_buckets = vegaBuckets(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
    cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);

disp('--- Vega-bucket sensitivity of the Certificate ---')
disp('Date | Vega (Notional = 100) | Vega (EUR) |')
for i = 1:length(vega_dates)
    disp([datestr(vega_dates(i), formatData), ' | ', num2str(vega_buckets(i) * 100), ' | ', ...
        num2str(vega_buckets(i) * Notional)]);
end
disp('--- --- ---')
 
%% Plot the vega bucket sensitivities

% plot_vega_buckets(vega_buckets, ttms);
% plot_vega_buckets(vega_buckets * Notional, ttms);

%% Coarse grained buckets

%% Delta

% define the buckets
buckets = [2,5,10,15];

% compute the coarse grained buckets for the delta
coarse_delta_buckets = deltaCoarseBuckets(dates(1), delta_dates, delta_buckets);

% Print the results of the coarse grained delta buckets
disp('--- Coarse grained delta buckets for the Certificate ---')
disp('Bucket (years) | DV01 (Notional = 100) | DV01 (EUR) |')
for i = 1:length(buckets)
    disp([num2str(buckets(i)), ' | ', num2str(coarse_delta_buckets(i) * 100), ' | ', ...
        num2str(coarse_delta_buckets(i) * Notional)]);
end
disp('--- --- ---')

%% Plot the coarse grained delta buckets

% plot_coarse_delta_buckets(buckets, coarse_delta_buckets*Notional);

%% Compute the sensitivities for the swaps for hedging

% compute the delta for the 4 swaps (2, 5, 10, 15 years)
swapRates = zeros(length(buckets), 1);

% initialize the matrix for the coarse grained delta buckets of the swaps
coarse_delta_buckets_swaps = zeros(length(buckets));
for i = 1:length(buckets)

    % compute the swap rate
    swapRates(i) = swapPricer(0, buckets(i), discounts, dates);

    % compute the bucket sensitivity
    [delta_buckets_swap, delta_dates_swap] = ...
        deltaBucketsSwap(datesSet, ratesSet, quoted_dates, swapRates(i), buckets(i));

    % compute the coarse grained buckets for the delta and store them
    coarse_delta_buckets_swaps(:,i) = deltaCoarseBuckets(dates(1), delta_dates_swap, delta_buckets_swap);
    
end

%% Plot the coarse grained delta buckets for the swaps

% plot_coarse_delta_buckets_swaps(buckets, coarse_delta_buckets_swaps*Notional);

%% Hedging of the Delta of the certificate

% compute the weights for the hedging with the swaps
delta_weights = HedgeCertificateDeltaBuckets(buckets, coarse_delta_buckets, coarse_delta_buckets_swaps, false);

disp('--- Swap notionals for Delta Hedging the Certificate by Maturity ---')
disp('Swap Maturity | Notional (EUR) |')
for i = 1:length(buckets)
    disp([num2str(buckets(i)), ' | ', num2str(delta_weights(i) * Notional)]);
end
disp('--- --- ---')

%% Plot the weights of the hedging

% plot_hedging_weights(buckets, delta_weights*Notional, 'Notional of the Delta hedging');

%% Vega hedging of certificate with ATM 5y Cap

% Hedge the Vega with an ATM 5y Cap (strike = ATM 5y Swap rate same conventions), and hedge the total portfolio 

% find the ATM 5y Cap strike
strike_5y = swapPricer(0, 5, discounts, dates) * 100;

% compute the vega for the ATM 5y Cap
[vega_5y_cap, cap_price_5y] = totalVegaCap(strike_5y, 5, spot_vols, spot_ttms, mkt_vols, ttms, strikes, discounts, dates);

% completely hedge the vega of the certificate with the ATM 5y Cap
weight_5y_cap = - total_vega / vega_5y_cap;

disp('--- Vega hedging with 5y ATM Cap ---')
disp(['Notional: ', num2str(weight_5y_cap * Notional), ' EUR']);
disp('--- --- ---')

%% Compute the coarse grained delta buckets of the 5y Cap

% compute the delta for the 5y Cap
[delta_dates_5y_cap, delta_buckets_5y_cap] = ...
    deltaBucketsCap(cap_price_5y, datesSet, ratesSet, quoted_dates, strike_5y, ...
    5, spot_vols, spot_ttms, strikes);

% compute the coarse grained buckets for the delta
coarse_delta_buckets_5y_cap = deltaCoarseBuckets(dates(1), delta_dates_5y_cap, delta_buckets_5y_cap);

%% New portfolio delta hedging

% compute the new portfolio delta
portfolio_delta = coarse_delta_buckets + weight_5y_cap * coarse_delta_buckets_5y_cap;

% compute the new weights for the hedging with the swaps
delta_weights_with_cap = HedgeCertificateDeltaBuckets(buckets, portfolio_delta, coarse_delta_buckets_swaps, false);

%% Print the results of the delta hedging and vega hedging

disp('--- Swap notionals for Delta Hedging the Certificate and 5y Cap by Maturity ---')
disp('Swap Maturity | Notional (EUR) |')
for i = 1:length(buckets)
    disp([num2str(buckets(i)), ' | ', num2str(delta_weights_with_cap(i) * Notional)]);
end
disp('--- --- ---')

%% Coarse grained vega buckets

% define the buckets
buckets_years = [5, 15];

% compute the coarse grained buckets for the vega of the certificate
coarse_vega_buckets = vegaCoarseBuckets(dates(1), vega_dates, vega_buckets);

%% Print the results of the vega buckets

disp('--- Coarse grained vega buckets for the Certificate ---')
disp('Bucket (years) | Notional (EUR) |')
for i = 1:length(buckets_years)
    disp([num2str(buckets_years(i)), ' | ', num2str(coarse_vega_buckets(i) * Notional)]);
end
disp('--- --- ---')

%% Coarse grained vega buckets for the 5y Cap

% compute the vega buckets for 5y cap
vega_buckets_5y_cap = vegaBucketsCap(cap_price_5y, strike_5y, 5, mkt_vols, ttms, strikes, discounts, dates);

% compute the coarse grained buckets for the vega of the 5y cap
coarse_vega_buckets_5y_cap = vegaCoarseBuckets(dates(1), vega_dates, vega_buckets_5y_cap);

%% Coarse grained vega buckets for the 15y Cap

% find the ATM 15y Cap strike
strike_15y = swapPricer(0, 15, discounts, dates) * 100;

% compute the price of the 15y Cap
[~, cap_price_15y] = totalVegaCap(strike_15y, 15, spot_vols, spot_ttms, mkt_vols, ttms, strikes, discounts, dates);

% compute the vega for the ATM 15y Cap
vega_buckets_15y_cap = vegaBucketsCap(cap_price_15y, strike_15y, 15, mkt_vols, ttms, strikes, discounts, dates);

% compute the coarse grained buckets for the vega of the 15y cap
coarse_vega_buckets_15y_cap = vegaCoarseBuckets(dates(1), vega_dates, vega_buckets_15y_cap);

%% Hedge the vega of the certificate with the 5y Cap and the 15y Cap

vega_weights = zeros(2, 1);

% compute the portfolio vega
portfolio_vega = coarse_vega_buckets;

% find the weight to perfectly hedge the bucket using the corresponding cap
vega_weights(2) = - portfolio_vega(2) / coarse_vega_buckets_15y_cap(2);

% update the portfolio vega
portfolio_vega = portfolio_vega + vega_weights(2) * coarse_vega_buckets_15y_cap;

% find the weight to perfectly hedge the bucket using the corresponding cap
vega_weights(1) = - portfolio_vega(1) / coarse_vega_buckets_5y_cap(1);

% update the portfolio vega
portfolio_vega = portfolio_vega + vega_weights(1) * coarse_vega_buckets_5y_cap;

disp('--- Notionals for Coarse Buckets Vega hedging with 5y and 15y ATM Cap ---')
disp(['Notional 5y Cap: ', num2str(vega_weights(1) * Notional)]);
disp(['Notional 15y Cap: ', num2str(vega_weights(2) * Notional)]);

toc;