%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Contest 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace
clc;
clear all;
close all;

tic

%% filename in .csv
% DATA

% IG Corporate Bonds DEsk





% % Compute Z_scores
% %
% % INPUT
% % zeroCurve : ZC bond data [yearfrac ; rates]
% % couponSchedule: cash flows bond [yearfrac; cash flow]
% % dirtyPrice = market dirty price
% %
% % OUTPUT
% % z = z_score

% ZC_curve = ;
% couponSchedule = ;
% dirtyPrice = ;

% z = zscore(ZC_curve, couponSchedule, dirtyPrice)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % price in one year of the risky bond
% %
% % INPUTS:
% %   IG_cf_schedule_2y: cash flow schedule of the IG bond
% %   Q: transition matrix
% %   ZC_curve: zero coupon curve
% %   Recovery: recovery rate
% %
% % OUTPUTS:
% %   FV: fair value of the bond in one year for IG, HY and Default

% IG_cf_schedule_2y = ;
% Q = ;
% ZC_curve = ;
% Recovery = ;

% FV = FV_risky_bond(IG_cf_schedule_2y, Q, ZC_curve, Recovery) 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Correlation factor for the IRB model
% %
% % INPUTS
% % PD: probability of default
% %
% % OUTPUTS
% % R: correlation factor

%PD = ;

% R= R_IRB(PD)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % QMATRIX  Derive the market-implied rating transition matrix based on the
% % available market data (hazard rates)
% %
% % INPUT
% % IG_h: hazard rates for the investment grade bonds
% % HY_h: hazard rates for the high yield bonds
% %
% % OUTPUT
% % Q: market-implied rating transition matrix

% IG_h = ;
% HY_h = ;

% Q = Qmatrix(IG_h,HY_h)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute the present value of a risky bond using the zero curve
% and the Z-score
% 
% INPUTS
% zScore: the Z-score of the bond
% couponSchedule: a matrix with the coupon dates and values
% ZC_curve: a matrix with the zero curve
% 
% OUTPUT
%   PV : Dirty price for a given risky bond (from scalar Z-spread)
% 
% z = 30 * 10^-4;
% cf_schedule = [0.5 0; 1.0 105.00] ;
% ZC_curve = [0.25 0.064; 0.5 0.063; 2.0 0.042];
% 
% PV_1y = PV_risky_bond_Z(z, cf_schedule, ZC_curve)
% 
% z = 50 * 10^-4;
% cf_schedule = [0.5 0; 1.0 5.5; 1.5 0; 2.0 105.50];
% 
% PV_2Y = PV_risky_bond_Z(z, cf_schedule, ZC_curve)
% 
% PV = (PV_1y * 4 + PV_2Y * 7 ) / 100
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PV_risky_bond_h: computes the present value of a risky bond using the
% survival probabilities (derived from the hazard rates) and the zero
% coupon curve
%
% INPUT
%   cf_schedule: a 2xN matrix with the cash flow schedule of the bond
%   h_curve: a 2xN matrix with the hazard rates
%   ZC_curve: a 2xN matrix with the zero coupon curve
%   R: recovery rate
%
% OUTPUT
%   PV : Dirty price for a given risky bond (from hazard rate curve - fixed recoovery rate R)

ZC_curve = [0.25 0.064; 0.5 0.063; 2.0 0.042];

h_curve = [1 , 59*10^-4; 2, 679*10^-4];
R = 0.25;

% One year bond IG
Notional_1y = 4*10^6;
cf_schedule = [0.5 0; 1.0 105.00];
IG_Bond_dirty_price_1y = 98.5234;
 
PV_1 = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R)

% Two year bond IG
Notional_2y = 7 * 10^6;
cf_schedule = [0.5 0; 1.0 5.5; 1.5 0; 2.0 105.50];
IG_Bond_dirty_price_2y = 98.5336;

PV_2 = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R)

PV = (PV_1 * 4 + PV_2 * 7 ) / 100

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Bootstrap the hazard 
% % 
% % INPUT
% % zeroCurve : ZC bond data [yearfrac ; rates]
% % R : recovery rate
% % dirtyPrice_1y : corporate market price 1year bond
% % dirtyPrice_2y : corporate market price 2year bond
% % couponSchedule_1y : cash flows 1year bond [yearfrac; cash flow]
% % couponSchedule_2y : cash flows 2year bond [yearfrac; cash flow]
% %
% % OUTPUT
% % h_curve : Hazard curve [yearfrac; h]

% zeroCurve = ;
% R = ;
% dirtyPrice_1y = ;
% dirtyPrice_2y = ;
% couponSchedule_1y = ;
% couponSchedule_2y = ;


% h_curve = hazardCurve(zeroCurve, R, dirtyPrice_1y, dirtyPrice_2y, ...
%     couponSchedule_1y, couponSchedule_2y)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc