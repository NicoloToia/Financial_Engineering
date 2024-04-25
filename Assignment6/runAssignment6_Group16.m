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

% plot_rates_discounts(dates, discounts, zeroRates);

%% Obtain the Cap Prices from the market data via Bachelier formula

% Cap prices
mkt_cap_prices = MarketCapPrices(ttms, strikes, mkt_vols, discounts, dates)

%% Compute the spot volatilites

spot_vols = spotVols(mkt_cap_prices, ttms, strikes, mkt_vols, discounts, dates);

% compute the exercise ttms
exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12*ttms(end)-3)';
ACT_360 = 2;
spot_ttms = yearfrac(dates(1), exercise_dates, ACT_360);

%% Plot the spot volatilities surface against the flat volatilities

plot_vols(spot_vols, spot_ttms, mkt_vols, ttms, strikes);

%% Price the certificate

spol_A = 2;
fixed_rate_B = 3;
spol_B = 1.1;
cap_5y = 4.3;
cap_10y = 4.6;
cap_15y = 5.1;

X = computeUpfront(spot_vols, spot_ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, discounts, dates)

toc;