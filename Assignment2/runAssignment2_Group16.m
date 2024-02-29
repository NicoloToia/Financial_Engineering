% runAssignment2
%  group X, AY20ZZ-20ZZ
% Computes Euribor 3m bootstrap with a single-curve model
%
% This is just code structure: it should be completed & modified (TBM)
%
% to run:
% > runAssignment2_TBM

clear;
close all;
clc;
% set the clock to find the time of execution
tic;
%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data
% This fuction works on Windows OS. Pay attention on other OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);
%% Bootstrap
% dates includes SettlementDate as first date
[dates, discounts]=bootstrap(datesSet, ratesSet);

readDates = datetime(dates, 'ConvertFrom', 'datenum');

%% Compute Zero Rates

zeroRates = zeroRates(dates, discounts);

%% Plot Results

plotresult(dates, discounts, zeroRates);


%% Point 2
% S is 2.8173%
fixedRate = 2.8173/100;

[ratesSet_shift] = shift_rates(ratesSet);

[dates, discounts_DV01]=bootstrap(datesSet, ratesSet_shift);

% target is 6 years
setDate = datesSet.settlement;
fixedLegPaymentDates = datesSet.swaps(1:6);

[DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, discounts,discounts_DV01)

MacD = sensCouponBond(setDate, fixedLegPaymentDates, fixedRate, dates, discounts)

%% Point 3

% price the IB coupon bond 6y with coupon rate equal to the mid-market swap at 7y
date_6Y = datesSet.swaps(6);
date_7Y = datesSet.swaps(7);
EU_30_360 = 6;
DF6Y = discounts(dates==date_6Y);
DF7Y = discounts(dates==date_7Y);
delta = yearfrac(date_6Y, date_7Y, EU_30_360);
S_7Y = 0.5 * (ratesSet.swaps(7,1) + ratesSet.swaps(7,2));

IB_couponBond = 1 + DF6Y - (1 + S_7Y * delta) * DF7Y

% direct calculation


toc;