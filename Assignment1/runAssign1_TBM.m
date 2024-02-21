% Assignment_1
%  Group 16, AA2023-2024
%

%% Clear the workspace
clear
close all
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

optionPrice = zeros(3,1);

for i=1:3
    pricingMode = i; % 1 ClosedFormula, 2 CRR, 3 Monte Carlo
    optionPrice(i) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
end

fprintf(['\nOPTION PRICE \n' ...
        'BlackPrice :   %.4f \n'],optionPrice(1));
fprintf('CRRPrice   :   %.4f \n',optionPrice(2));
fprintf('MCPrice    :   %.4f \n',optionPrice(3));

%% Point b

% spread is 1 bp
spread = 10^-4;
% find optimal M for CRR
optionPriceBlack = EuropeanOptionClosed(F0,K,B,TTM,sigma,flag);
<<<<<<< Updated upstream
M_CRR = findMCRR (optionPriceBlack, F0, K, B, TTM, sigma, flag, spread)

% find optimal M for MC
[M_optimal]=OptimalM_MC(F0,K,B,TTM,sigma, spread)
=======
M_CRR = findMCRR (optionPriceBlack, F0, K, B, TTM, sigma, flag, spread);

fprintf(['\nOPTIMAL M FOR CRR \n' ...
        'Number of intervals CRR :   %.d \n'],M_CRR);

>>>>>>> Stashed changes
%% Point c

% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 

% Plot the results of CRR
figure
loglog(nCRR,errCRR)
title('CRR')
hold on
loglog(nCRR, 1./nCRR)
% cutoff
loglog(nCRR, spread * ones(length(nCRR),1))

% Plot the results of MC
figure
loglog(nMC,stdEstim)
title('MC')
hold on 
loglog(nMC,1./sqrt(nMC))
% cutoff
loglog(nMC, spread * ones(length(nMC),1))

%% Point d

% set the barrier
KI = 1.3;

% store the prices
optionPriceKI = zeros(3,1);

% % closed formula
% optionPriceKI(1) = EuropeanOptionKIClosed(F0,K,KI,B,TTM,sigma);
% % CRR
% optionPriceKI(2) = EuropeanOptionKICRR(F0,K,KI,B,TTM,sigma, 1000);
% % monte carlo
% optionPriceKI(3) = EuropeanOptionKIMC(F0,K,KI,B,TTM,sigma,1000000);

% closed formula
optionPriceKI(1) = EuropeanOptionKIClosed(F0,K,KI,B,TTM,sigma);
% CRR
optionPriceKI(2) = EuropeanOptionKICRR(F0,K,KI,B,TTM,sigma, M);
% monte carlo
optionPriceKI(3) = EuropeanOptionKIMC(F0,K,KI,B,TTM,sigma,M);

fprintf(['\nOPTION PRICE KI \n' ...
        'ClosedPriceKI  :   %.4f \n'],optionPriceKI(1));
fprintf('CRRPriceKI     :   %.4f \n',optionPriceKI(2));
fprintf('MCPriceKI      :   %.4f \n',optionPriceKI(3));

%% KI Option

%% KI Option Gamma

%% Antithetic Variables
