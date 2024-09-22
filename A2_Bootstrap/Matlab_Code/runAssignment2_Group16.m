% runAssignment2_Group16
%  group 16, AY2023-2024
% Computes Euribor 3m bootstrap with a single-curve model
%
% to run:
% > runAssignment2_Group16

% clear workspace
clear all;
close all;
clc;

% set the clock to find the time of execution
tic;

%% Settings

addpath Data\
addpath Functions\
formatData='dd/mm/yyyy'; % Pay attention to your computer settings 

%% Read market data

% This function works on Windows OS. Pay attention on other OS.
[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);

% Faster loading dataSet
% load datesSet.mat
% load ratesSet.mat

%% Bootstrap
tic;
% Bootstrap the discount factors from the market data
% dates includes SettlementDate as first date (for more detail see function bootstrap)
[dates, discounts] = bootstrap(datesSet, ratesSet);

% save the results
save('discounts.mat', 'dates', 'discounts');

% Covert dates in a more readable format
readDates = datetime(dates, 'ConvertFrom', 'datenum');

%% Compute Zero Rates

% Compute zero rates from the discounts founded in the bootstrap
zeroRates = zeroRates(dates, discounts);

%% Plot Results

% Plot the zero rates curve and discounts curve obtained by bootstrap technique
plotresult(dates, discounts, zeroRates);


%% Sensitivities (Point 2)
% With the discount curve obtained above compute the sensitivities for a portfolio composed 
% only by one single swap, a 6y plain vanilla IR swap vs Euribor 3m 
% with a fixed rate S = 2.8173%

fixedRate = 2.8173/100;

% Parallel shift of market rates by 1bp
[ratesSet_shift] = shift_rates(ratesSet);

% Boostrap the curve with the shifted rates
[dates, discounts_DV01] = bootstrap(datesSet, ratesSet_shift);

% Set dates, target date is 6 years (i.e. end calculation)
setDate = datesSet.settlement;
fixedLegPaymentDates = datesSet.swaps(1:6);

% Compute sensitivities: DV01, BPV, DV01_z (see function sensSwap for more detail)
[DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01);

% Compute Macaulay duration
MacD = sensCouponBond(setDate, fixedLegPaymentDates, fixedRate, dates, discounts);

% Display results
fprintf(['\nSWAP SENSITIVITIES \n' ...
        'DV01   :   %.4e \n'],DV01);
fprintf('BPV    :   %.4e \n',BPV);
fprintf('DV01_z :   %.4e \n',DV01_z); 
fprintf('MacD   :   %.4f \n',MacD)

%% Point 3
% Price the IB coupon bond 6y with coupon rate equal to the mid-market swap at 7y

% Set dates
date_6Y = datesSet.swaps(6);    % end of calculation for 6Y IB coupon Bond
date_7Y = datesSet.swaps(7);    % end of calculation for swap 7Y

% yearfrac convention
EU_30_360 = 6;
delta = yearfrac(date_6Y, date_7Y, EU_30_360);  % delta(t6,t7)

% Set discounts
DF6Y = discounts(dates==date_6Y);   % B(t0,t6)
DF7Y = discounts(dates==date_7Y);   % B(t0,t7)
% Find mid-market swap rate 7Y
S_7Y = 0.5 * (ratesSet.swaps(7,1) + ratesSet.swaps(7,2));

% Compute the price of the bond by shortcut (see the report for more details)
IB_couponBond = 1 + DF6Y - (1 + S_7Y * delta) * DF7Y;

% Display results
fprintf(['\nPRICE IB COUPON BOND 6Y WITH RATE EQUAL TO MID SWAP 7Y \n' ...
        'IB Coupon Bond Price  :   %.4f \n\n'],IB_couponBond);

%%
toc;
