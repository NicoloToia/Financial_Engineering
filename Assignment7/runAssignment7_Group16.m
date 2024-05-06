% runAssignment7_Group16
%  group 16, AY2023-2024
% 
%
% to run:
% > runAssignment7_Group16

clc
clear all
close all

% fix random seed
rng(42) 

%% Load data for exercise 1
principal = 1e8; % 100 million
settlement_date = datenum(datetime(2008, 2, 19)); % 19-Feb-2008
maturity = 2; % 2 years
EU_30_360 = 6; % 30/360 day count convention   
coupon_1y = 0.06; 
coupon_2y = 0.02; 
strike = 3200; 
trigger = 0.06;

%coupon dates each year
couponDates = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(1:2)';
couponDates(~isbusday(couponDates, eurCalendar())) = busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);

%discount factors
load('discounts.mat');

% calibrated parameters for alpha = 1/2
Sigma = 0.10493;
Kappa = 1.3089;
Eta = 12.7545;




