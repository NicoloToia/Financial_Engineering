% Assignment_1
%  Group 16, AA2023-2024
% commento nick

%% Clear the workspace
clear all
close all
warning('off','all')
clc

%% Fix the random seed
rng(42);

%% Pricing parameters
S0=1;
K=1;
r=0.03;
TTM=1/4; 
sigma=0.22;
flag=1; % flag:  1 call, -1 put
d=0.06;

%% Quantity of interest
B=exp(-r*TTM); % Discount

%% Pricing 
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model

%% Point a
M=100; % M = simulations for MC, steps for CRR;

OptionPrice = zeros(3,1);

for i=1:3
    pricingMode = i; % 1 ClosedFormula, 2 CRR, 3 Monte Carlo
    OptionPrice(i) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
end

%% Point b

% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 

%% Point c

% Plot the results of CRR
figure
title('CRR')
loglog(nCRR,errCRR)
hold on
loglog(nCRR,1/nCRR)

% Plot the results of MC
figure
title('MC')
loglog(M,stdEstim)
hold on
loglog(M,1./sqrt(M))

%% KI Option

%% KI Option Gamma

%% Antithetic Variables