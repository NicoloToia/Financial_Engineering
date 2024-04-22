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
[ttms, strikes, Vols] = readVolData('Caps_vol_20-2-24');

%% Construct the swap dates
 
swapDates = datetime(datesSet.settlement, 'ConvertFrom', 'datenum') + calyears(1:50)';
swapDates(~isbusday(swapDates, eurCalendar())) = busdate(swapDates(~isbusday(swapDates, eurCalendar())), 'modifiedfollow', eurCalendar());
swapDates = datenum(swapDates);

%% Interpolate the mid market rates for swaps from the market data

% find the actually quoted swap dates
datesSet.swaps = [swapDates(1:12); swapDates(15:5:30); swapDates(40:10:50)];

% Interpolation of the mid market rates for swaps
ratesSet.swaps = interp1(datesSet.swaps, ratesSet.swaps, swapDates, 'spline', 'extrap');
datesSet.swaps = swapDates;

%% Bootstrap the discount factors from the market data

% Bootstrap the discount factors from the market data
[dates, discounts] = bootstrap(datesSet, ratesSet);

% compute zero rates from the discounts founded in the bootstrap
zeroRates = zeroRates(dates, discounts);

%% Plot Results

plotresult(dates, discounts, zeroRates);

%% Obtain the Cap Prices from the market data via Bachelier formula

% Cap prices
capPrices = capPrices(ttms, strikes, Vols, ratesSet, swapDates);

toc;