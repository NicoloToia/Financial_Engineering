% runAssignment2
%  group X, AY20ZZ-20ZZ
% Computes Euribor 3m bootstrap with a single-curve model
%
% This is just code structure: it should be completed & modified (TBM)
%
% to run:
% > runAssignment2_TBM

clear all;
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

[dates, discounts, zeroRates]=bootstrap(datesSet, ratesSet);

% plot
figure 
yyaxis left
plot(dates(2:end), discounts(2:end), 'r-s')
ylabel('Discount Factors')
yyaxis right
plot(dates(2:end), 100*zeroRates(2:end), 'b-o')
ylabel('Zero Rates')
xlabel('Date')
title('Zero Rates and Discount Factors')
legend('Discount Factors', 'Zero Rates')
grid on

%% Compute Zero Rates
% TBM

%% Plot Results

%discount curve
% TBM

%zero-rates
% TBM
toc;