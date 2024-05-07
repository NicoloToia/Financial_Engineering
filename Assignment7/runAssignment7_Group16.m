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

principal = 100e6; % 100 million
ttm = 2; % 2 years
coupon_1y = 6 / 100;
coupon_2y = 2 / 100;
s_A = 1.3 / 100;
strike = 3200; 
trigger = 6 / 100;

% initial price and parameters for the NIG
S_0 = cSelect.reference;
d = cSelect.dividend;
sigma = calibrated_parameters(1);
kappa = calibrated_parameters(2);
eta = calibrated_parameters(3);

%% Compute the upfront payment

X = computeUpfront(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, discounts, dates);

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate ---')
disp(['The upfront payment is: ', num2str(X/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X), ' EUR']);
disp('--- --- ---')

%% Tree

% compute the reset dates for the tree
tree_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(1:10)';
% convert to business date
tree_dates(~isbusday(tree_dates,eurCalendar())) = ...
    busdate(tree_dates(~isbusday(tree_dates,eurCalendar())),'modifiedfollow',eurCalendar());
tree_dates = datenum(tree_dates);

% data for the tree

a = 0.11;
sigma = 0.008;
dt = 1/200;
ttm = 10;

[l_max, mu, trinomial_tree] = buildTrinomialTree(a, sigma, dt, ttm);

% print the trinomial tree
disp('--- Trinomial Tree ---')
disp('The trinomial tree is:')
for i = 1:ttm/dt
    disp(['At time ', num2str(i), ' the tree is: ', num2str(trinomial_tree{i}')]);
end

% find wh