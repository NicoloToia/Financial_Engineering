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
principal = 100*1e9; % 100 billion
settlement_date = datenum(datetime(2008, 2, 19)); % 19-Feb-2008
maturity = 2; % 2 years
ACT_30_360 = 2; % 30/360 day count convention   
coupon_1y = 0.06; 
coupon_2y = 0.02; 
strike = 3200; 
trigger = 0.06;

%coupon dates each year
couponDates = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(1:2)';
couponDates(~isbusday(couponDates, eurCalendar())) = busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);



