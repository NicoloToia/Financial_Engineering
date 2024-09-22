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
rng(42) % the answer to everything in the universe

%% Add the directories to the path

addpath('Data');
addpath('Discounts');
addpath('FFT');
addpath('pricing')
addpath('Tree')
addpath('schemes')
addpath('jamshidian')

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

%% Compute the the upfront via Closed Integral method

X = Compute_Upfront_Closed(S_0, d, strike, 2, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, ...
    discounts, dates, alpha);

% print the upfront
disp('--- Upfront payment of the Certificate computed via close integral formula ---')
disp(['The upfront payment is: ', num2str(X/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X), ' EUR']);
disp('--- --- ---')

%% Compute upfront MC and NIG

N = 1e6;

X = computeUpfrontMCV(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, ...
    discounts, dates, alpha, N);

% print the upfront payment
disp('--- Upfront payment of the Certificate MCV ---')
disp(['The upfront payment is: ', num2str(X/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X), ' EUR']);
disp('--- --- ---')

%% Compute the upfront payment by NIG and FFT

X_NIG = computeUpfrontFFT(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, ...
     discounts, dates, alpha);

% print the upfront payment percentage
disp('--- Upfront payment of the Certificate ---')
disp(['The upfront payment is: ', num2str(X_NIG/principal*100), '%']);
disp(['The upfront payment is: ', num2str(X_NIG), ' EUR']);
disp('--- --- ---')

%% Compute the upfront payment via Variance Gamma

% uncomment to perform the VG method

% % load('calibrated_parameters_gamma.mat')
% % 
% % alpha = 0;
% % sigma = calibrated_parameters_gamma(1);
% % kappa = calibrated_parameters_gamma(2);
% % eta = calibrated_parameters_gamma(3);
% % 
% % X_VG = computeUpfrontFFT(S_0, d, strike, ttm, principal, coupon_1y, coupon_2y, s_A, sigma, kappa, eta, discounts, dates, alpha);
% % 
% % % print the upfront payment percentage
% % disp('--- Upfront payment of the Certificate VG---')
% % disp(['The upfront payment is: ', num2str(X_VG/principal*100), '%']);
% % disp(['The upfront payment is: ', num2str(X_VG), ' EUR']);
% % disp('--- --- ---')


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

%% 3y bond

ttm = 3; 
N = 1e6; 
coupons = [0.06, 0.06, 0.02];
alpha_IC = 5/100; % confidence level at 5%

[X_3y, IC_3y] = price3y(S_0, d, strike, ttm, alpha, sigma, kappa, eta, s_A, N, discounts, dates, ...
    principal, coupons, alpha_IC);
disp ('--- Pricing bond with expiry 3y ---')
disp(['The upfront payment is: ', num2str(X_3y/principal*100), '%']);
disp(['The 95% confidence interval is [', num2str(IC_3y(1)/principal*100), '%, ' , ...
    num2str(IC_3y(2)/principal*100), '%]'])
disp(['The price of the 3y bond is   : ', num2str(X_3y)]);
disp(['The 95% confidence interval is [', num2str(IC_3y(1)), ', ' , ...
    num2str(IC_3y(2)), '] EUR'])
disp(['IC width is: ', num2str(IC_3y(2) - IC_3y(1)), ' EUR'])
disp('--- --- ---')

%% BERMUDAN SWAPTION BY TRINOMIAL TREE

% import the data
strike = 5/100; 
ttm = 10; % 10 years time to maturity
% the option can be exercided every 1 year from the 2 year (non-call 2)
% HULL-WHITE parameters
a = 11/100;
sigma = 0.8/100;

% compute the reset dates for the tree from 2 to 9 years
reset_dates = (datenum(dates(1)) + (1:ttm)*365)';


% set the number of steps in each interval
N_steps_in_1y = 70;
% N_steps_in_1y = 4;
N_steps = N_steps_in_1y * ttm;

% find the interval length dt
dt = 1/N_steps_in_1y;

% find the dates in each node
node_dates = (datenum(dates(1)) + dt*(0:N_steps)*365)';
% find the integer for the dates
node_dates = round(node_dates);


% compute sigma hat and mu hat
mu_hat = 1- exp(-a*dt);
sigma_hat = sigma * sqrt((1-exp(-2*a*dt))/(2*a));
% find the jump D_x in the tree
D_x = sigma_hat * sqrt(3);
% find l_max and l_min (symmetric)
l_max = ceil((1-sqrt(2/3))/mu_hat);
l_min = -l_max;
% vector of position in the tree (l)     
l = l_max:-1:l_min;

% build the vector of x_i based on the indicator l
x = l * D_x;

% compute the forward discount factor in each reset date
% fwdDE_reset = zeros(2*l_max+1, length(resetDates));

discounts_reset = intExtDF(discounts, dates, reset_dates);
BPV = zeros(2*l_max+1, length(reset_dates)-1);
swaps_rate = zeros(2*l_max+1, length(reset_dates)-1);
intr_value = zeros(2*l_max+1, length(reset_dates)-1);
intr_value1 = zeros(2*l_max+1, length(reset_dates)-1);

for i = 2:9

    fwdDF_ttm = discounts_reset(end)/discounts_reset(i);
    fwdDF_present = fwdDF_ttm*exp( -x * (sigmaHJM(a, sigma, (ttm-i), 0)/sigma) - 0.5 * IntHJM(a, sigma, i, 0, (ttm-i) ) );
    float_leg = (1 - fwdDF_present)';
    % find the BPV
    fwdDF_dt = discounts_reset(i+1:end)/discounts_reset(i);
    BPV = zeros(2*l_max+1, 1);

    for j = 1:ttm-i
        fwdDF_present_dt = fwdDF_dt(j) * exp( -x * (sigmaHJM(a, sigma, j, 0)/sigma) - 0.5 * IntHJM(a, sigma, i, 0, j));
        BPV = BPV + fwdDF_present_dt'; % is an yearly bpv so delta_i is always 1
    end

%     BPV = BPV(:,i)
    % find the swap rates
    swaps_rate(:,i) = float_leg ./ BPV;

    % intr_value1(:,i) = max(0, 1 - BPV * strike + fwdDF_present');
    intr_value(:,i) = BPV .* max(0, swaps_rate(:,i) - strike);

end

%% Build the tree with the stochastic dicount in each node for the step (i, i+1)

discounts_node = intExtDF(discounts, dates, node_dates);
discounts_node(1) = 1;
fwd_discount_nodes = discounts_node(2:end)./discounts_node(1:end-1);
sigma_star = (sigma/a) * sqrt(dt - 2 *( (1 - exp(-a*dt)) / a ) + (1 - exp(-2*a*dt)) / (2*a) );

% initialize
fwdDF_present = zeros(2*l_max+1, N_steps);

for i = 1:N_steps

    fwdDF_present(:,i) = fwd_discount_nodes(i)*exp( -x * (sigmaHJM(a, sigma, dt, 0)/sigma) - 0.5 * IntHJM(a, sigma, i, 0, dt) );
    
end

value = zeros(2*l_max + 1, (ttm)*N_steps_in_1y+1);

for i = (ttm)*N_steps_in_1y:-1:1

    % compute the continuation value
    for j = 1:2*l_max+1

        if j == 1
            value(j,i) = C_contvalue(l_max, mu_hat, sigma, sigma_star, a, dt, fwdDF_present(:,i), value(:,i+1), D_x, x);

        elseif j == 2*l_max+1
            value(j,i) = B_contvalue(l_min, mu_hat, sigma, sigma_star, a, dt, fwdDF_present(:,i), value(:,i+1), D_x, x);

        else
            value(j,i) = A_contvalue(l_max - j + 1, j, mu_hat, sigma, sigma_star, a, dt, fwdDF_present(:,i), value(:,i+1), D_x, x);            
        end

    end
    
    if find(reset_dates == node_dates(i))
    
        % find the index of the reset date
        index = find(reset_dates == node_dates(i));
    
        continuation_value = max(intr_value(:,index), value(:,i));
    
        value(:,i) = continuation_value;
    
    else
        continue
    end

    
    % value(:,i) = continuation_value;

end

% display the value of the swaption
disp('--- Bermudan swaption value ---')
disp(['The value of the Bermudan swaption is: ', num2str(value(l_max+1,1))]);

%% Jamshidian Formula

% expiry of the swaption and start
alfa = 9;
ttm = 10;
% stike of the swaption
strike_swaption = 5 / 100;

% swaption from 2 to 10
Jamshidian(alfa, ttm, strike_swaption, a, sigma, discounts, dates);

% display the value of the swaption
disp('--- Swaption value ---')
disp(['The value of the swaption (9y-10y) is: ', num2str(Jamshidian(alfa, ttm, strike_swaption, a, sigma, discounts, dates))]);
disp('--- --- ---')

% find the price for each alpha
alphas = 2:ttm-1;

% initialize the vector of prices
prices = zeros(length(alphas), 1);

for i = 1:length(alphas)

    prices(i) = Jamshidian(alphas(i), ttm, strike_swaption, a, sigma, discounts, dates);

end

% find upper and lower bound
upper_bound = sum(prices);
lower_bound = min(prices);