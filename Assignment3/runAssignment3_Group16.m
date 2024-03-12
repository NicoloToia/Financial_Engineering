% runAssignment3_Group16
%  group 16, AY2023-2024
% Compute the 
%
% to run:
% > runAssignment2_Group16

% clear workspace
clear all;
close all;
clc;

format long

% set the clock to find the time of execution
tic;

%% Settings

formatData='dd/MM/yyyy'; % Pay attention to your computer settings 

%% Load the data from the previous assignment
% discounts defines both dates and the discount factors
load('discounts.mat');

% compute the floating leg dates
settlementDate = datetime(dates(1), 'ConvertFrom', 'datenum');
swap3yDate = datetime(dates(13), 'ConvertFrom', 'datenum');
floatDates = datenum(settlementDate+calmonths(3):calmonths(3):swap3yDate);
% get the business days with the modified following convention
% we specify 0 to avoid American holidays, indeed 21/02/2011 was Presidents' day (NYSE is closed)
floatDates(~isbusday(floatDates)) = busdate(floatDates(~isbusday(floatDates)), "modifiedfollow", 0);

% compute the discount factors at the float dates
discountsFloat = intExtDF(discounts, dates, floatDates);

% find the fixed leg dates
swap1yDate = datenum(datetime('19/02/2009', 'InputFormat', formatData));
% fixed leg dates
fixedDates = [swap1yDate; dates(12); dates(13)];
% fixed leg discount factors
discountsFixed = intExtDF(discounts, dates, fixedDates);

% coupon parameters
C_bar_0 = 101/100;
C_bar = 3.9/100;

% compute the price of the corresponding IB coupon bond
ACT_360 = 2;
deltas = yearfrac([dates(1); fixedDates(1:end-1)], fixedDates, ACT_360);
C0 = C_bar * deltas' * discountsFixed + discountsFixed(end);

% compute the bpv using the floating leg
ACT_360 = 2;
deltas = yearfrac([dates(1); floatDates(1:end-1)], floatDates, ACT_360);
BPV_float = deltas' * discountsFloat;

% compute the asset swap spread (expressed in basis points)
S_asw = (C0 - C_bar_0) / BPV_float;
S_asw = S_asw * 10000;
% display the result
disp(['The asset swap spread is: ', num2str(S_asw), ' basis points'])

%% Point 2

% recovery rate
R = 0.4;

% take the swap dates 1y -> 5y, 7y
datesCDS = [swap1yDate; dates(12:15); dates(17)];
% spreads in basis points
spreadsCDS = [ 29, 32, 35, 39, 40, 41] / 10000;


%% Point 2.a

% create the spline for the complete set of dates
completeDates = [swap1yDate; dates(12:17)];
% use cubic spline to interpolate the spreads
spreadsCDS = spline(datesCDS, spreadsCDS, completeDates);

% plot the spreads 
figure
plot(completeDates, spreadsCDS, 'o-')
title('CDS Spreads')
xlabel('Dates')
ylabel('Spreads')


%% Point 2.b
t0 = dates(1);
% Andrea
% [datesCDS, survProbs, intensities] =  bootstrapCDS(dates, discounts, [t0 ; completeDates], spreadsCDS, 1, R)
% Jacopo
[datesCDS, survProbs, intensities] =  bootstrapCDS_v2(dates, discounts, completeDates, spreadsCDS, 1, R)

PlotIntensities(datesCDS, intensities,t0)

toc