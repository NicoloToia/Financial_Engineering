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

% Interpolation of the mid market rates for swaps using spline
ACT_365 = 3;
delta_swaps_set = yearfrac(datesSet.settlement, datesSet.swaps, ACT_365);
delta_swaps = yearfrac(datesSet.settlement, swapDates, ACT_365);
ratesSet.swaps = interp1(delta_swaps_set, ratesSet.swaps, delta_swaps, 'spline');
datesSet.swaps = swapDates;

%% Bootstrap the discount factors from the market data

% Bootstrap the discount factors from the market data
[dates, discounts] = bootstrap(datesSet, ratesSet);

% compute zero rates from the discounts founded in the bootstrap
zeroRates = zeroRates(dates, discounts);

%% Plot Results

% plotresult(dates, discounts, zeroRates);

%% Obtain the Cap Prices from the market data via Bachelier formula

% Cap prices
mkt_cap_prices = MarketCapPrices(ttms, strikes, mkt_vols, discounts, dates)

%% Compute the spot volatilites

spot_vols = spotVols(mkt_cap_prices, ttms, strikes, mkt_vols, discounts, dates)

% plot the spot volatilities surface against the flat volatilities
% spot vols in blue, market vols in red
surf(strikes, ttms, spot_vols, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'blue');
hold on
surf(strikes, ttms, mkt_vols, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'red');
xlabel('Strikes');
ylabel('TTMs');
zlabel('Volatilities');
title('Spot Volatilities vs Market Volatilities');
legend('Spot Volatilities', 'Market Volatilities');

toc;