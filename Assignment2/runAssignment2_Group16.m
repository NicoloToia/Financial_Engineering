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
tic;
%% Settings
formatData='dd/mm/yyyy'; %Pay attention to your computer settings 

%% Read market data
% This fuction works on Windows OS. Pay attention on other OS.

[datesSet, ratesSet] = readExcelData('MktData_CurveBootstrap', formatData);
%% Bootstrap
% dates includes SettlementDate as first date
[dates, discounts]=bootstrap(datesSet, ratesSet);


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

[DV01, BPV, DV01_z] = sensSwap(setDate, fixedLegPaymentDates, fixedRate, dates, ...
    discounts,discounts_DV01)

MacD = sensCouponBond(setDate, fixedLegPaymentDates, fixedRate, dates, discounts)

toc;