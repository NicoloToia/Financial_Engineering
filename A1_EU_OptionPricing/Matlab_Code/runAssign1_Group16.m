% Assignment1_Group16
% Group 16, AA2023-2024
%
%% Clear the workspace
clear; close all; clc
addpath Functions\

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
B=exp(-r*TTM); % Discount Factor

%% Pricing (Point a)
F0=S0*exp(-d*TTM)/B;     % Forward in G&C Model

%   Price the option, 
%   considering an underlying price equal to 1 Euro (i.e a derivative Notional of 1 Mln Euro): 
%   i) via blkprice Matlab function; 
%   ii) with a CRR tree approach; 
%   iii) with a Monte-Carlo (MC) approach. 

M=100; % M =  number of simulations for MC & time steps for CRR;

optionPrice = zeros(3,1); % Inizialize option price vector 

% Price by the function EuropeanOptionPrice
for i=1:3
    pricingMode = i; % 1 ClosedFormula, 2 CRR, 3 Monte Carlo
    optionPrice(i) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode,M,flag);
end

% Display results
fprintf(['\nOPTION PRICE \n' ...
        'BlackPrice :   %.4f \n'],optionPrice(1));
fprintf('CRRPrice   :   %.4f \n',optionPrice(2));
fprintf('MCPrice    :   %.4f \n',optionPrice(3));

% notional amount
notional = 1e6;
fprintf(['\nOPTION PRICE FOR 1 Mln CONTRACTS \n' ...
        'BlackPrice :   %.4f \n'],notional * optionPrice(1));
fprintf('CRRPrice   :   %.4f \n',notional * optionPrice(2));
fprintf('MCPrice    :   %.4f \n',notional * optionPrice(3));

%% Errors Rescaling (Point b & c)

% Consider M, as the number of intervals in CRR and as the number of 
% simulations in the MC. Focus on a call.Select M according to the criteria
% mentioned in the class.
% Show that the numerical errors for a call rescale with M as 1/𝑀 for CRR 
% and as 1/√𝑀 for MC. 

% Plot Errors for CRR varing number of steps
[nCRR,errCRR]=PlotErrorCRR(F0,K,B,TTM,sigma);
% Plot Errors for MC varing number of simulations N 
[nMC,stdEstim]=PlotErrorMC(F0,K,B,TTM,sigma); 

% Spread is 1 bp
spread = 10^-4;

% Find the optimal M for CRR (between the admissible choiche of vector M(1x10) ) 
M_CRR = nCRR(find(errCRR < spread,1));
% Find the optimal M for MC (between the admissible choiche of vector M(1x20) )
M_MC = nMC(find(stdEstim < spread,1));

% Find best estimate CRR
M_CRR_opt = findMCRR (optionPrice(1), F0, K, B, TTM, sigma, flag, spread);

% Display results
fprintf(['\nBEST ADMISSIBLE M FOR CRR \n' ...
        'Number of intervals CRR :   %.d \n'],M_CRR);
fprintf(['\nBEST M FOR MC \n' ...
        'Number of intervals CRR :   %.d \n'],M_MC);
fprintf('\nOPTIMAL M FOR CRR:   %d \n',M_CRR_opt);

%% KI Option (Point d)

%   Price a European Call Option with European barrier at €1.3 (up&in) and same parameters with 
%   the two numerical techniques (tree & MC). Does it exist a closed formula also in this case? If yes, 
%   compare the results. 

% Set the barrier
KI = 1.3;

% Set a suitable number of time steps for CRR and simulations for MC
M = 1000;  

% Inizialize option price KI vector
optionPriceKI = zeros(3,1);

% Compute the price by means of the two numerical methods and the closed
% formula
% Closed formula
optionPriceKI(1) = EuropeanOptionKIClosed(F0,K,KI,B,TTM,sigma);
% CRR
optionPriceKI(2) = EuropeanOptionKICRR(F0,K,KI,B,TTM,sigma,M);
% Monte Carlo
optionPriceKI(3) = EuropeanOptionKIMC(F0,K,KI,B,TTM,sigma,M);

% DIsplay the results
fprintf(['\nOPTION PRICE KI \n' ...
        'ClosedPriceKI  :   %.4f \n'],optionPriceKI(1));
fprintf('CRRPriceKI     :   %.4f \n',optionPriceKI(2));
fprintf('MCPriceKI      :   %.4f \n',optionPriceKI(3));

%% KI Option Vega (Point e)

%   For this barrier option, plot the Vega (possibly using both the closed formula and a numerical 
%   estimate) with the underlying price in the range 0.70 Euro and 1.5 Euro. Comment the results. 

% Set underlying price range
S_start = 0.7;
S_end = 1.5;
% Set a grid oft the underlying price
rangeS0 = linspace(S_start,S_end,100);

% Compute the corresponding forward prices
rangeF0 = rangeS0*exp(-d*TTM)/B;

% Set M the number of MC simulations & time step of CRR
M = 10000;

% Closed formula vegas
vegasClosed = zeros(length(rangeS0),1);

for i = 1:length(rangeS0)
    vegasClosed(i) = VegaKI(rangeF0(i),K,KI,B,TTM,sigma,M,3);
end

% MC vegas
vegasMCEstim = zeros(length(rangeS0),1);

for i = 1:length(rangeS0)
    vegasMCEstim(i) = VegaKI(rangeF0(i),K,KI,B,TTM,sigma,M_MC,2);
end

% CRR vegas
vegasCRR = zeros(length(rangeS0),1);

for i = 1:length(rangeS0)
    vegasCRR(i) = VegaKI(rangeF0(i),K,KI,B,TTM,sigma,M,1);
end

% Plot the results
figure
subplot(1,2,1)
plot(rangeS0,vegasClosed)
title('Vega Closed Formula')
xlabel('S0')
ylabel('Vega')
hold on
plot(rangeS0,vegasMCEstim)
title('Vega Monte Carlo')
xlabel('S0')
ylabel('Vega')

legend('Closed','MC')

subplot(1,2,2)
plot(rangeS0,vegasClosed)
title('Vega Closed Formula')
xlabel('S0')
ylabel('Vega')
hold on 
plot(rangeS0,vegasCRR)
title('Vega CRR')
xlabel('S0')
ylabel('Vega')

legend('Closed','CRR')

%% MC Error Estim

vegasMCEstim = zeros(length(rangeS0),1);
for i = 1:length(rangeS0)
    vegasMCEstim(i) = VegaMCestim(rangeF0(i),K,KI,B,TTM,sigma,M);
end

figure
plot(rangeS0,vegasMCEstim)
hold on
plot(rangeS0,vegasClosed)
title('Vega Monte Carlo')
xlabel('S0')
ylabel('Vega')

legend('Vega estimator','Closed')

%% Antithetic Variables (Point f)

% Does antithetic variables technique (Hull 2009, Ch.19.7) reduce MC error of point b.? 

% Plot Errors for MCAV varing number of simulations N 
[nMCAV,stdEstimAV]=PlotErrorMCAV(F0,K,B,TTM,sigma);


% Plot the results of MCAV
figure
loglog(nMCAV,stdEstimAV)
title('MC Antithetic Variable Error')
xlabel('M'); ylabel('errorMC AV')
hold on 
loglog(nMCAV,1./sqrt(nMCAV))
loglog(nMCAV, spread * ones(length(nMCAV),1))   % cutoff based on the spread
loglog(nMC,stdEstim)

legend('MC AV','1/sqrt(M)','cutoff','MC')

MAV=132000; %   optimal number of iterations for MCAV
%   price for the European call option vith the AC method
optionPriceAV = EuropeanOptionMCAV(F0,K,B,TTM,sigma,MAV,flag);
fprintf('\nOPTION PRICE AV : %.4f \n',optionPriceAV);

%% Bermudan Option (Point g)

%   Price also -with the Tree- a Bermudan option, where the holder has also the right to exercise 
%   the option at the end of every month, obtaining the stock at the strike price. 
M = 1000;

[OptionPriceBermudan,~] = BermudanOptionCRR(F0, K, B, TTM, sigma, d, M);

% Display results
fprintf('\nCRRPriceBermudan      :   %.4f \n',OptionPriceBermudan);

%%  Vary the Dividend Yield (Point h)

%   Pricing the Bermudan option, vary the dividend yield between 0% and 6% and compare 
%   with the corresponding European price. Discuss the results.

% Set range of dividend yields between 0% and 6%
d_range=0:0.005:0.06;

OptionPriceBermudan = zeros(length(d_range),1);
OptionPriceAM = zeros(length(d_range),1);
OptionCRR = zeros(length(d_range),1);

for i=1:length(d_range)
    % recompute the forward price
    F0=S0*exp(-d_range(i)*TTM)/B;
    [OptionPriceBermudan(i), OptionPriceAM(i)] = BermudanOptionCRR(F0, K, B, TTM, sigma, d_range(i), M);
    OptionCRR(i) = EuropeanOptionCRR(F0, K, B, TTM, sigma, M, flag);
end

% Plot the results
figure
plot(d_range, OptionPriceBermudan, '-xb', 'LineWidth', 1); 
hold on
plot(d_range, OptionPriceAM, '-xg', 'LineWidth', 1); 
hold on; 
plot(d_range, OptionCRR, '-or', 'LineWidth', 1);
xlabel('Dividend Yield'); 
ylabel('Option Price'); 
title('Comparison of Bermudan Option Price and European Option Price'); 
legend('Bermudan Option Price', 'Pseudo-American Option', 'European Option Price'); 
