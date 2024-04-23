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

% % compute the first cap price for the first strike and first year

% payment_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12)';
% % move to business days if needed
% payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
%     busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar())
% payment_dates = datenum(payment_dates);

% datestr(dates(1))

% deltas = yearfrac([dates(1);payment_dates(1:end-1)], payment_dates, 2)

% % compute the discount factors
% DF = intExtDF(discounts, dates, payment_dates)

% % forward discounts and forward Libor
% fwdDiscounts = DF ./ [1;intExtDF(discounts, dates, payment_dates(1:end-1))]
% fwdLibor = 1 ./ deltas .* (1 ./ fwdDiscounts - 1)

% sigma = Vols(1,1)
% K = strikes(1)
% ttm = ttms(1)

% % skip the first caplet

% % price the second caplet (6m)
% delta_ttm = yearfrac(dates(1), payment_dates(1), 2)
% d_n = (fwdLibor(2) - (1+K/100)*fwdLibor(2)) / (sigma * sqrt(delta_ttm));

% Caplet_6m = DF(2) * deltas(2) * ( (fwdLibor(2) - (1+K/100) * fwdLibor(2)) * normcdf(d_n) + sigma * sqrt(delta_ttm) * normpdf(d_n))

% CapletBachelier(fwdLibor(2), (1+K/100)*fwdLibor(2), sigma, deltas(2), payment_dates(2), dates(1), DF(2))

% Cap prices
mkt_cap_prices = MarketCapPrices(ttms, strikes, Vols, discounts, dates);

toc;