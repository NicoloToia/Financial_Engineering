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
addpath('CapPricing');
addpath('Sensitivities');
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
    cap_5y, cap_10y, cap_15y, discounts, dates);

% print the upfront payment percentage
disp(['The upfront payment is: ', num2str(X*100), '%']);

%% Delta-bucket sensitivity

% compute the delta-bucket sensitivity
[delta_dates, delta_buckets] = deltaBuckets(datesSet, ratesSet, dates, spot_vols, spot_ttms, strikes, X, dates(1), spol_A, ...
    fixed_rate_B, spol_B, cap_5y, cap_10y, cap_15y);

%% Plot the delta-bucket sensitivities

plot_delta_buckets(delta_dates, delta_buckets);

%% Total Vega

total_vega = total_Vega(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, discounts, dates);

% print the total vega
disp(['The total vega is: ', num2str(total_vega), ' bp']);

%% Vega bucket sensitivity

% if vegas files is present, load it, otherwise compute it
if isfile('Data/vega_buckets.mat')
    load('vega_buckets.mat');
    % retransform the vegas
    vega_buckets = vega_buckets / 10^6;
    % divide by the shift
    vega_buckets = vega_buckets / 0.0001;
else
    vega_buckets = vegaBuckets(mkt_vols, ttms, strikes, X, spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    save('vega_buckets.mat', 'vega_buckets');
end

% print the vega buckets sensitivities as a table

% create the row names from the ttms
row_names = cell(1, length(ttms));
for i = 1:length(ttms)
    row_names{i} = ['T = ', num2str(ttms(i)), 'y'];
end

% create the column names from the strikes
column_names = cell(1, length(strikes));
for i = 1:length(strikes)
    column_names{i} = ['K = ', num2str(strikes(i))];
end

% create the table
vega_buckets_table = array2table(vega_buckets, 'RowNames', row_names, 'VariableNames', column_names);

% print the table
disp(vega_buckets_table);

%% Plot the vega bucket sensitivities

plot_vega_buckets(vega_buckets, ttms, strikes);

toc;