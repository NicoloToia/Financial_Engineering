% Assignment_1
%  Group X, AA2021-2022
%
%  TBM (To Be Modified): Modify & Add where needed

%% Pricing parameters
S0=1;
K=1;
r=0.04;
TTM=1/6; 
sigma=0.30;
flag=1; % flag:  1 call, -1 put
d=0.05;

%% Quantity of interest
B=exp(-r*TTM); % Discount

%% Pricing 
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model

%TBM: Modify with a cicle
pricingMode = 1; % 1 ClosedFormula, 2 CRR, 3 Monte Carlo
M=100; % M = simulations for MC, steps for CRR;
OptionPrice = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag)

%% Errors Rescaling 

% plot Errors for CRR varing number of steps
% Note: both functions plot also the Errors of interest as side-effect 
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);

% plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 

%% KI Option

%% KI Option Gamma

%% Antithetic Variables