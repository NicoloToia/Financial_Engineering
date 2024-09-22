%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case RM 2
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT data are stored as:
%
%   1. ZC_curve: table of zero-coupon rates (continuous compounding)
%       Column #1: maturity (year frac)
%       Column #2: MID rate
% 
%   2. XX_cf_schedule_Yy: table of cash flows of corp. bonds with 
%       XX rating and Yy maturity
%       Column #1: cash flow date(year frac)
%       Column #2: cash flow amount (US $)
%      alsongside with they market prices (dirty)
% 
%   legend: XX -> IG = Investment Grade
%           XX -> HY = High Yield
%
%    OUTPUT data are stored as:
%
%   3. IG_FV: Column vector with 1-year forward value of IG bond,
%      each row correaponding to a rating 
%      (row #1: IG; row #2: HY; row #3: Defaulted)
%      Example: FV = [102.00; 98.00; 45.00]
%   5. A set of scalar variable with self-explanatory names
%   
%   Required functions' template:
%   1. Dirty 1y fwd. prices for a given risky bond 
%      output: column vector (see above specs.)
%     fwd. prices are derived from transition matrix, fixed recovery rate R
%     function [ FV ] = FV_risky_bond(cf_schedule, Q, ZC_curve, R)  
%   2. IRB correlation (large corporate-sov.)for exposure with given PD 
%   function [R] = R_IRB(PD) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fix seed 
rng(2);
clc;
clear all;
addpath Functions\
% format("bank")
format long

% Zero Coupon Curve
ZC_curve = [0.25 0.054; 0.5 0.053; 2.0 0.0487];
 
% Two year bond IG
IG_cf_schedule_2y = [0.5 2.7407; 1.0 2.7407; 1.5 2.7407; 2.0 102.7407];
IG_Bond_dirty_price_2y = 100.00000;

% Rating transition matrix
Q = [0.7790	0.2160	0.0050;
     0.4079	0.5529	0.0392;
     0.0000	0.0000	1.0000];

%Recovery rate (\pi according to Schonbucher)
Recovery = 0.40;                        

% Number of issuers in portfolio
N_issuers = 200;    
N_issuers = 20; % Uncomment to assess concentration risk

% AVR correlation (Schonbucher 10.12) - uncomment for discussion
R=0.00;
% R=0.12;
% R=0.24; 
% R = R_IRB(Q(1,3));
rho = sqrt(R);
% disp('––– Part II Q3: Basel II correlation function –––')
% fprintf('IG correlation: %.3f \n', rho)
% disp(' ')


% Minimum number of Monte Carlo Scenarios
N = 100000;
% N = N * 100;       % Uncomment for convergence test (10,000,000 scenarios)

%% Q1: CreditMetrics' FV (let's start with all bonds rated IG)
FV = FV_risky_bond(IG_cf_schedule_2y, Q, ZC_curve, Recovery); 
E_FV = Q(1,1)*FV(1) + Q(1,2)*FV(2) + Q(1,3)*FV(3);

disp('––– Part I Q1: Present Value in a years’time of the IG bond –––')
fprintf('1y fwds: %.2f IG; %.2f HY; %.2f Def \n',FV )
fprintf('Present Value in a years time: %.2f\n',E_FV)
disp(' ')


%% Barriers and simulated P/L for a single IG issuer
bD = norminv(Q(1,3));                       % Barrier to Def.
bd = norminv(Q(1,2)+Q(1,3));                % Barrier to down (HY)
bu = 1.0e6;                                     % Barrier to up (never)
L_D = (FV(3) - E_FV)/N_issuers;                 % Loss if default
L_d = (FV(2) - E_FV)/N_issuers;                 % Loss if down
L_i = (FV(1) - E_FV)/N_issuers;                 % Profit if unchanged
L_u = 0;                                        % Profit if up 

%% Monte Carlo
% Generate N standard 1-d normal variables (macro-factor)
Y = randn(N,1);                      % Realization of the single factor

%% Idiosynchratic Monte Carlo
% zeros(N,1);        % double: Realization of the AVR 
scen_D = zeros(N,1);   % int.: Count of default events in each MC scenario
scen_d = zeros(N,1);   % int.: Count of "down" events in each MC scenario
scen_i = zeros(N,1);   % int.: Count of "identical rating" events in each MC scenario
scen_u = zeros(N,1);   % int.: Count of "up" events in each MC scenario
h = waitbar(0, 'Loading...'); % Initialize waitbar
for i = 1:N_issuers
    V = rho * Y + sqrt(1 - rho^2) * randn(N,1); % Realization of the idiosyncratic factor
    % Update the event counts
    scen_D = scen_D + (V <= bD);
    scen_d = scen_d + (V <= bd & V > bD);
    scen_i = scen_i + (V > bd & V <= bu);
    scen_u = scen_u + (V > bu);
    % Update waitbar
    waitbar(i/N_issuers, h, sprintf('Loading... %d/%d', i, N_issuers));
end

close(h); % Close waitbar after loop completes


%% Test on MC accuracy

% test on marginal distributions
D_count = sum(scen_D)/N;        % Average number of default events in portfolio
d_count = sum(scen_d)/N;        % Average number of "down" events in portfolio
i_count = sum(scen_i)/N;        % Average number of "identical rating" events in portfolio
u_count = sum(scen_u)/N;        % Average number of "up" events in portfolio
% Are MC results in agreement with the input rating transition matrix?
pD_simulated =D_count / N_issuers;
pd_simulated =d_count / N_issuers;
pi_simulated =i_count / N_issuers;
pu_simulated =u_count / N_issuers;
% ...if simulated probabilities are not satisfactory, maybe N is too low

fprintf('––– Test of the accuracy of the MC simulation (%.0f names)–––\n', N_issuers)
fprintf('PD         (simulated) %.4f; PD         (input) %.4f\n',[pD_simulated,Q(1,3)])
fprintf('down prob. (simulated) %.4f; down prob. (input) %.4f \n',[pd_simulated,Q(1,2)])
fprintf('stay prob. (simulated) %.4f; stay prob. (input) %.4f \n',[pi_simulated,Q(1,1)])
fprintf('up prob.   (simulated) %.4f; up prob.   (input) %.4f \n',[pu_simulated,0])
disp(' ')

%% Questions 2 and 3 (and discussion): Credit VaR

% Default only case (DO)
EL = L_D * pD_simulated + L_d * pd_simulated + L_i * pi_simulated + L_u * pu_simulated;
L_y = quantile(-L_D * scen_D,0.999);
VaR_DO = L_y - EL;

% Default and migration case (DM)
L_y = quantile(-L_D * scen_D - L_d * scen_d - L_i * scen_i,0.999);
VaR_DM = L_y - EL;

fprintf('––– Part I Q2/3: Credit VaR with AVR correlation %.3f (%.0f names)–––\n', [R,N_issuers])
fprintf('VaR - default only %.2f \n',VaR_DO)
fprintf('VaR - default and migration %.2f \n',VaR_DM)
disp(' ')
