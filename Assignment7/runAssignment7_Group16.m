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
addpath('pricing')
addpath('Tree')
addpath('schemes')


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
alpha = 1/2;
sigma = calibrated_parameters(1);
kappa = calibrated_parameters(2);
eta = calibrated_parameters(3);

%% Compute upfront mc

N = 1e6;

X = computeUpfrontMCV(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, ...
    discounts, dates, alpha, N);

% print the upfront payment
disp('--- Upfront payment of the Certificate ---')
disp(['The upfront payment is: ', num2str(X/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X), ' EUR']);
disp('--- --- ---')

%%
X = Compute_Upfront_Closed(S_0, d, strike, 2, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, ...
    discounts, dates, alpha);

% print the upfront
disp('--- Upfront payment of the Certificate computed via close integral formula ---')
disp(['The upfront payment is: ', num2str(X/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X), ' EUR']);
disp('--- --- ---')



%% Compute the upfront payment

X_NIG = computeUpfrontFFT(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, discounts, dates, alpha);

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate ---')
disp(['The upfront payment is: ', num2str(X_NIG/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X_NIG), ' EUR']);
disp('--- --- ---')

%% Compute the upfront payment via Variance Gamma

% % % load('calibrated_parameters_gamma.mat')
% % % 
% % % alpha = 0;
% % % sigma = calibrated_parameters_gamma(1);
% % % kappa = calibrated_parameters_gamma(2);
% % % eta = calibrated_parameters_gamma(3);
% % % 
% % % X_VG = computeUpfrontFFT(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, discounts, dates, alpha);
% % % 
% % % % print the upfront payment percentage
% % % disp('--- Upfront payment of the Certificate VG---')
% % % disp(['The upfront payment is: ', num2str(X_VG/principal*100), '%']);
% % % disp(['The upfront payment is: ', num2str(X_VG), ' EUR']);
% % % disp('--- --- ---')


%% Black con Skew

% save the data to use
quoted_strikes = cSelect.strikes;
quoted_vols = cSelect.surface;

flag = 2; % Black with skew

X_black_skew = computeUpfrontSkew(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, discounts, dates, ...
    quoted_strikes, quoted_vols, flag);

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate via Black adj skew---')
disp(['The upfront payment is: ', num2str(X_black_skew/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X_black_skew), ' EUR']);
disp('--- --- ---')



%% Black price without skew (no digital risk)

flag = 1; % Black without skew

X_black = computeUpfrontSkew(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, discounts, dates, ...
    quoted_strikes, quoted_vols, flag);

% find the error in the price wrt the one with skew and the one with NIG
black_vs_blackSkew = abs(X_black - X_black_skew) / X_black_skew * 100;
black_vs_NIG = abs(X_black - X_NIG) / X_NIG * 100;
blackSkew_vs_NIG = abs(X_black_skew - X_NIG) / X_NIG * 100;

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate via Black---')
disp(['The upfront payment is: ', num2str(X_black/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X_black), ' EUR']);
disp('--- --- ---')

disp('--- The Error using Black model (%) ---')
disp(['The error between Black without skew and adj Black   : ', num2str(black_vs_blackSkew)]);
disp(['The error between Black without skew and NIG         : ', num2str(black_vs_NIG)]);
disp(['The error between Black with skew and NIG            : ', num2str(blackSkew_vs_NIG)]);
disp('--- --- ---')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% disp('--- Trinomial Tree ---')
% disp('The trinomial tree is:')
% for i = 1:ttm/dt
%     disp(['At time ', num2str(i), ' the tree is: ', num2str(trinomial_tree{i}')]);
% end

% find wh