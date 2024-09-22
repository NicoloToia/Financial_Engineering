%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case RM 1
% Group 16
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
%      alongside with they market prices (dirty)
% 
%   legend: XX -> IG = Investment Grade
%           XX -> HY = High Yield
% 
%    OUTPUT data are stored as:
%
%   1. XX_h_curve: Table of piece-wise hazard rates
%      Maturities are year fractions. No maturity 0 is accepted:
%      Example: IG_h_curve = [1.00 0.01; 2.00 0.02; 5.00 0.05]
%      Hazard rate is constant at 1% from t=0 to t=1y, then it sharply 
%      increases to 2% from t=1y to t=2y, then sharlpy increases to 5%
%      and is constant till 5y
%   2. Q(3,3): 1y rating transition matrix
%   3. A set of scalar variable with self-explanatory names
%   
%   Required functions' template:
%   1. Dirty price for a given risky bond (from hazard rate curve - fixed recoovery rate R)
%   function [ PV ] = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R) 
%   2. Dirty price for a given risky bond (from scalar Z-spread)
%   function [ PV ] = PV_risky_bond_Z( z, cf_schedule, ZC_curve)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
addpath Functions\
tic
format("bank")

% Zero Coupon Curve
ZC_curve = [0.25 0.054; 0.5 0.053; 2.0 0.046];

% One year bond IG
IG_cf_schedule_1y = [0.5 2.5; 1.0 102.500];
IG_Bond_dirty_price_1y = 99.5705;
 
% Two year bond IG
IG_cf_schedule_2y = [0.5 2.75; 1.0 2.75; 1.5 2.75; 2.0 102.75];
IG_Bond_dirty_price_2y = 100.5420;

% One year bond HY
HY_cf_schedule_1y = [0.5 2.250; 1.0 102.250];
HY_Bond_dirty_price_1y = 97.0445;

% Two year bond HY
HY_cf_schedule_2y = [0.5 2.375; 1.0 2.375; 1.5 2.375; 2.0 102.375];
HY_Bond_dirty_price_2y = 96.4160;

R = 0.4;                        %Recovery rate (\pi according to Schonbucher)

%%
% Part I Q1: Bootstrap of hazard rate curves (two, each one corresponding to a rating class) 

% Bootstrap the hazard curve
IG_h_curve = hazardCurve(ZC_curve, R, IG_Bond_dirty_price_1y, ...
    IG_Bond_dirty_price_2y, IG_cf_schedule_1y, IG_cf_schedule_2y);
    
HY_h_curve = hazardCurve(ZC_curve, R, HY_Bond_dirty_price_1y, ...
    HY_Bond_dirty_price_2y, HY_cf_schedule_1y, HY_cf_schedule_2y);

% Display hazard curve
disp('––– Part I Q1: Hazard rate piecewise curves –––')
fprintf('IG Hazard rate (bp): %.1f 1y; %.1f 2y \n',IG_h_curve(:,2)*10000 )
fprintf('HY Hazard rate (bp): %.1f 1y; %.1f 2y \n',HY_h_curve(:,2)*10000 )
disp(' ')

% Check: Model price of the bonds should match the market prices
IG_Bond_dirty_price_1y_check = PV_risky_bond_h(IG_cf_schedule_1y, IG_h_curve, ZC_curve,R);
IG_Bond_dirty_price_2y_check = PV_risky_bond_h(IG_cf_schedule_2y, IG_h_curve, ZC_curve,R);
HY_Bond_dirty_price_1y_check = PV_risky_bond_h(HY_cf_schedule_1y, HY_h_curve, ZC_curve,R);
HY_Bond_dirty_price_2y_check = PV_risky_bond_h(HY_cf_schedule_2y, HY_h_curve, ZC_curve,R);

% Display model prices vs market prices
disp('––– Check whether model prices (hazard rate) match input market prices –––')
fprintf('IG 1y model price: %.4f (Mkt: %.4f)\n' , IG_Bond_dirty_price_1y_check,IG_Bond_dirty_price_1y) 
fprintf('IG 2y model price: %.4f (Mkt: %.4f)\n', IG_Bond_dirty_price_2y_check,IG_Bond_dirty_price_2y)
fprintf('HY 1y model price: %.4f (Mkt: %.4f)\n', HY_Bond_dirty_price_1y_check,HY_Bond_dirty_price_1y)
fprintf('HY 2y model price: %.4f (Mkt: %.4f)\n', HY_Bond_dirty_price_2y_check,HY_Bond_dirty_price_2y)
disp(' ')

%%
% Part I Q2: Z-spread (four scalars, each one corresponding to a bond)

% Compute the zero score
IG_z_1y = zscore(ZC_curve, IG_cf_schedule_1y, IG_Bond_dirty_price_1y);
IG_z_2y = zscore(ZC_curve, IG_cf_schedule_2y, IG_Bond_dirty_price_2y);
HY_z_1y = zscore(ZC_curve, HY_cf_schedule_1y, HY_Bond_dirty_price_1y);
HY_z_2y = zscore(ZC_curve, HY_cf_schedule_2y, HY_Bond_dirty_price_2y);

IG_z_curve = [IG_z_1y; IG_z_2y];
HY_z_curve = [HY_z_1y; HY_z_2y];

% Display the zero scores results
disp('––– Part I Q2: Z-spread (scalar): basis points –––')
fprintf('IG 1y Z-spread: %.1f bp \n', IG_z_1y*10000) 
fprintf('IG 2y Z-spread: %.1f bp \n', IG_z_2y*10000)
fprintf('HY 1y Z-spread: %.1f bp \n', HY_z_1y*10000)
fprintf('HY 2y Z-spread: %.1f bp \n', HY_z_2y*10000)
disp(' ')

% Check: Model price of the bonds should match the market prices
IG_Bond_dirty_price_1y_check = PV_risky_bond_Z(IG_z_1y, IG_cf_schedule_1y, ZC_curve);
IG_Bond_dirty_price_2y_check = PV_risky_bond_Z(IG_z_2y, IG_cf_schedule_2y, ZC_curve);
HY_Bond_dirty_price_1y_check = PV_risky_bond_Z(HY_z_1y, HY_cf_schedule_1y, ZC_curve);
HY_Bond_dirty_price_2y_check = PV_risky_bond_Z(HY_z_2y, HY_cf_schedule_2y, ZC_curve);

% Display model prices vs market prices
disp('––– Check whether model prices (Z-spread) match input market prices –––')
fprintf('IG 1y model price: %.4f (Mkt: %.4f)\n' , IG_Bond_dirty_price_1y_check,IG_Bond_dirty_price_1y) 
fprintf('IG 2y model price: %.4f (Mkt: %.4f)\n', IG_Bond_dirty_price_2y_check,IG_Bond_dirty_price_2y)
fprintf('HY 1y model price: %.4f (Mkt: %.4f)\n', HY_Bond_dirty_price_1y_check,HY_Bond_dirty_price_1y)
fprintf('HY 2y model price: %.4f (Mkt: %.4f)\n', HY_Bond_dirty_price_2y_check,HY_Bond_dirty_price_2y)
disp(' ')

%% Part II Q1: Transition Matrix

% Compute the transition matrix
[Q] = Qmatrix(IG_h_curve, HY_h_curve);
% Print the transition matrix as a table
disp(array2table(Q, 'VariableNames', {'IG', 'HY', 'Default'}, 'RowNames', {'IG', 'HY', 'Default'}))

toc