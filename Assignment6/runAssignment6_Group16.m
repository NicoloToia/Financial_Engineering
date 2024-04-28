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
disp(['The upfront payment is: ', num2str(X*100), '%']);
Notional = 50 * 10^6;
disp(['The upfront payment is: ', num2str(X*Notional), ' EUR']);

%% Delta-bucket sensitivity

% compute the delta-bucket sensitivity
[delta_dates, delta_buckets] = deltaBuckets(datesSet, ratesSet, quoted_dates, spot_vols, spot_ttms, strikes, X, dates(1), spol_A, ...
    fixed_rate_B, spol_B, cap_rate_5y, cap_rate_10y, cap_rate_15y);

%% Plot the delta-bucket sensitivities

% plot_delta_buckets(delta_dates, delta_buckets);
% plot_delta_buckets(delta_dates, delta_buckets * Notional * 10^(-4));
% these values are expressed in bp, to obtain EUR values, multiply by the notional and by 1 bp

%% Total Vega

total_vega = total_Vega(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
    cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);
% print the total vega
disp(['The total vega is: ', num2str(total_vega), ' bp']);
disp(['The total vega is: ', num2str(total_vega * Notional * 10^(-4)), ' EUR']);

%% Vega bucket sensitivity

if ~isfile('Data/vega_buckets_vector.mat')
    % compute the vega bucket sensitivities
    vega_buckets = vegaBuckets(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
        cap_rate_5y, cap_rate_10y, cap_rate_15y, discounts, dates);
    
    save('vega_buckets_vector.mat', 'vega_buckets');
else
    load('vega_buckets_vector.mat');
end

%% Plot the vega bucket sensitivities

% plot_vega_buckets(vega_buckets, ttms);
% plot_vega_buckets(vega_buckets * Notional * 10^(-4), ttms);

%% Coarse grained buckets

%% Delta

buckets = [2,5,10,15];

% compute the coarse grained buckets for the delta
coarse_delta_buckets = deltaCoarseBuckets(dates(1), delta_dates, delta_buckets);

%% Plot the coarse grained delta buckets

% plot_coarse_delta_buckets(buckets, coarse_delta_buckets);

%% Compute the sensitivities for the swaps for hedging

% compute the delta for the 4 swaps (2, 5, 10, 15 years)
swapRates = zeros(length(buckets), 1);
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

% plot_coarse_delta_buckets_swaps(buckets, coarse_delta_buckets_swaps);

%% Hedging of the Delta of the certificate

delta_weights = HedgeCertificateDeltaBuckets(buckets, coarse_delta_buckets, coarse_delta_buckets_swaps, true);

%% Plot the weights of the hedging

% plot_hedging_weights(buckets, delta_weights, 'Weights of the Delta hedging');

%% Vega

% add the initial date to the vega dates
vega_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears([0;ttms]);
vega_dates(~isbusday(vega_dates, eurCalendar())) = ...
    busdate(vega_dates(~isbusday(vega_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
vega_dates = datenum(vega_dates);

% compute the coarse grained buckets for the vega (add zero to have the point of reference)
coarse_vega_buckets = vegaCoarseBuckets(vega_dates, [0;vega_buckets]);

%% Plot the coarse grained vega buckets

% plot_coarse_vega_buckets(buckets, coarse_vega_buckets);


%% Vega hedging of certificate with ATM 5y Cap

% Hedge the Vega with an ATM 5y Cap (strike = ATM 5y Swap rate same conventions), and hedge the total portfolio 

% find the ATM 5y Cap strike
strike_5y = swapPricer(0, 5, discounts, dates);

% compute the vega for the ATM 5y Cap
vega_5y_cap = totalVegaCap(strike_5y, 5, spot_vols, spot_ttms, mkt_vols, ttms, strikes, discounts, dates)

% completely hedge the vega of the certificate with the ATM 5y Cap
vega_weight_5y_cap = -total_vega / vega_5y_cap;

% print the vega weight
disp(['The vega weight for the 5y ATM Cap is: ', num2str(vega_weight_5y_cap)]);

%% Coarse grained vega buckets for the swaps

% aggregate the vega buckets for into 0-5y, 5y-15y

portfolio_coarse_buckets = [coarse_vega_buckets(1) + coarse_vega_buckets(2); 
    coarse_vega_buckets(3) + coarse_vega_buckets(4)];

% compute the vega of the 15 year cap 

% find the ATM 15y Cap strike
strike_15y = swapPricer(0, 15, discounts, dates);

% compute the vega for the ATM 15y Cap
% TODO: tecnically we should compute the bucket vega for the 15y cap
vega_15y_cap = vegaCap(strike_15y, 15, spot_vols, spot_ttms, mkt_vols, ttms, strikes, discounts, dates)

% start hedging the vega with the larger bucket
vega_weights = zeros(length(portfolio_coarse_buckets), 1);

% find the weight to perfectly hedge the bucket using the corresponding cap
vega_weights(2) = - portfolio_coarse_buckets(2) / vega_15y_cap;

% update the portfolio vega
portfolio_coarse_buckets = portfolio_coarse_buckets + vega_weights(2) * [0;vega_15y_cap];

% find the weight to perfectly hedge the bucket using the corresponding cap
vega_weights(1) = - portfolio_coarse_buckets(1) / vega_5y_cap;

% update the portfolio vega
portfolio_coarse_buckets = portfolio_coarse_buckets + vega_weights(1) * [0;vega_5y_cap];

toc;