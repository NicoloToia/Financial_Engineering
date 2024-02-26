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
% plot discount factors as green filled triangles
plot(dates(2:end), discounts(2:end), 'g-^', 'MarkerFaceColor', 'g')
ylabel('Discount Factors')
% make the y-axis go from 0 to 1
ylim([0 1])
% ticks every 0.2
yticks(0:0.2:1)
yyaxis right
% plot zero rates as blue filled diamonds (in percent)
plot(dates(2:end), zeroRates(2:end)*100, 'b-d', 'MarkerFaceColor', 'b')
% make the y-axis go from 2.5 to 5.0
ylim([2.5 5.0])
yticks(2.5:0.5:5.0)
ylabel('Zero Rates')
xlabel('Date')
title('Zero Rates and Discount Factors')
legend('Discount Factors', 'Zero Rates')
grid on
% x-axis in date format Month name and last two digits of year
datetick('x', 'mmm-yy')
% cut off dates and rates after the Jan-2042
% make axes black
ax = gca;
ax.XColor = 'k';
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

%% Compute Zero Rates
% TBM

%% Plot Results

%discount curve
% TBM

%zero-rates
% TBM
toc;