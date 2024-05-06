% runAssignment7_Group16
%  group 16, AY2023-2024
% 
%
% to run:
% > runAssignment7_Group16

clc
clear all
close all

tic

% fix random seed
rng(42) 

%% Add the directories to the path

addpath('Data');
addpath('Discounts');
addpath('FFT');

%% Load data for bootstrap and NIG model

%discount factors
load('discounts.mat');
load('calibrated_parameters.mat');
load('EuroStoxx_data.mat');

%% Data for exercise 1

principal = 1e8; % 100 million
settlement_date = dates(1); % 19 Feb 2008 
ttm = 2; % 2 years
EU_30_360 = 6; % 30/360 day count convention   
ACT_365 = 3;
coupon_1y = 0.06; 
coupon_2y = 0.02; 
strike = 3200; 
trigger = 0.06;

%coupon dates each year
couponDates = datetime(settlement_date, 'ConvertFrom', 'datenum') + calyears(1:ttm)';
couponDates(~isbusday(couponDates, eurCalendar())) = ...
    busdate(couponDates(~isbusday(couponDates, eurCalendar())), 'modifiedfollow', eurCalendar());
couponDates = datenum(couponDates);

% compute the discount factors at the payment dates
coupon_DF = intExtDF(discounts, dates, couponDates);

% find the date up to which to simulate the underlying
fixing_date = datetime(couponDates(1), 'ConvertFrom', 'datenum') - caldays(2);
fixing_DF = intExtDF(discounts, dates, datenum(fixing_date));
t = yearfrac(settlement_date, fixing_date, ACT_365);

% Compute the price of the digital as a call option spread

% initial data
S_0 = cSelect.reference;
d = cSelect.dividend;
F_0 = S_0 / fixing_DF * exp(-d * t);
epsilon = 1;
% parameters for the FFT
alpha = 1/2;
sigma = calibrated_parameters(1);
kappa = calibrated_parameters(2);
eta = calibrated_parameters(3);

% compute the desired log moneyness
x_1 = log(F_0/strike);
x_2 = log(F_0/(strike + epsilon));

%% Point 3.a: FFT method, alpha = 1/2

% compute the call prices with the FFT method
M_FFT = 15;
dz = 0.0025;
flag = 'FFT';
call_1 = callIntegral(fixing_DF, F_0, alpha, sigma, kappa, eta, t, x_1, M_FFT, dz, flag);
call_2 = callIntegral(fixing_DF, F_0, alpha, sigma, kappa, eta, t, x_2, M_FFT, dz, flag);

% take the difference to compute the digital
digital_FFT = (call_2 - call_1) / epsilon;
