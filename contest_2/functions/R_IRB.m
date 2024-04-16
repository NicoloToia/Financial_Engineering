function [R]= R_IRB(PD)
% Correlation factor for the IRB model
%
% INPUTS
% PD: probability of default
%
% OUTPUTS
% R: correlation factor

R_min = 0.12;
R_max = 0.24;
K = 50;

R = R_min * (1-exp(-K*PD)) / (1-exp(-K)) + R_max * (exp(-K*PD)) / (1-exp(-K));

end