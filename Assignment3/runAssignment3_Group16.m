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

% convert dates to datetime
dates = datetime(dates, 'ConvertFrom', 'datenum');

settlementDate = dates(1);
swap3yDate = dates(13);

% compute the floating leg dates
floatDates = (settlementDate+calmonths(3):calmonths(3):swap3yDate)';
% use the modified following convention to find the next business day
floatDates = busdate(floatDates, 'modifiedfollow', -1);
% compute the discount factors at the float dates
discountsFloat = intExtDF(discounts, datenum(dates), datenum(floatDates));

% find the fixed leg dates
swap1yDate = datetime('19/02/2009', 'InputFormat', formatData);
% fixed leg dates
fixedDates = [swap1yDate; dates(12); dates(13)];
% fixed leg discount factors
discountsFixed = intExtDF(discounts, datenum(dates), datenum(fixedDates));

% coupon parameters
C_bar_0 = 101/100;
C_bar = 3.9/100;

% compute the price of the corresponding IB coupon bond
ACT_360 = 2;
deltas = yearfrac([settlementDate; fixedDates(1:end-1)], fixedDates, ACT_360);
C0 = C_bar * deltas' * discountsFixed + discountsFixed(end);

% compute the bpv using the floating leg
ACT_360 = 2;
deltas = yearfrac([settlementDate; floatDates(1:end-1)], floatDates, ACT_360);
BPV_float = deltas' * discountsFloat;

% compute the asset swap spread
S_asw = (C0 - C_bar_0) / BPV_float;


toc