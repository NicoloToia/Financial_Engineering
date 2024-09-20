% runAssignment3_Group16
%  group 16, AY2023-2024
% Compute the 
%
% to run:
% > runAssignment2_Group16

% clear workspace
clear all;
close all;
clc;

% NOTE : for long computational time the number of iterations are set by
% default at 1e4, use a larger number if more precision is required
% (line:147)

%% Settings

formatData ='dd/MM/yyyy'; % Pay attention to your computer settings 

rng(42);   % Fix the random number generator ("the answer to Life, the Universe, and Everything")

tic;       % Set the clock to find the time of execution

%% Asset Swap (Point 1)

% Load the data from the previous assignment
% Discounts defines both dates and the discount factors
load('discounts.mat');

% Compute the floating leg dates
settlementDate = datetime(dates(1), 'ConvertFrom', 'datenum');
swap3yDate = datetime(dates(13), 'ConvertFrom', 'datenum');
% From the settlement date add 3 months until 3 years are reached
floatDates = datenum(settlementDate+calmonths(3):calmonths(3):swap3yDate)';
% Get the business days with the modified following convention
% We specify 0 to avoid American holidays, indeed 21/02/2011 was Presidents' day (NYSE is closed)
floatDates(~isbusday(floatDates)) = busdate(floatDates(~isbusday(floatDates)), "modifiedfollow", 0);

% Compute the discount factors at the float dates by intExDF function,
% which interpolates if necessary
discountsFloat = intExtDF(discounts, dates, floatDates);

% Find the fixed leg dates
swap1yDate = datenum(datetime('19/02/2009', 'InputFormat', formatData)); % 1 year swap date from data
% Fixed leg dates
fixedDates = [swap1yDate; dates(12); dates(13)];
% Fixed leg discount factors
discountsFixed = intExtDF(discounts, dates, fixedDates);

% Coupon parameters from data
C_bar_0 = 101/100;
C_bar = 3.9/100;

% Compute the price of the corresponding IB coupon bond
EU_30_360 = 6;
deltas = yearfrac([dates(1); fixedDates(1:end-1)], fixedDates, EU_30_360);
C0 = C_bar * deltas' * discountsFixed + discountsFixed(end);

% Compute the BPV using the floating leg
ACT_360 = 2;
deltas = yearfrac([dates(1); floatDates(1:end-1)], floatDates, ACT_360);
BPV_float = deltas' * discountsFloat;

% Compute the asset swap spread (expressed in basis points)
S_asw = (C0 - C_bar_0) / BPV_float;
S_asw = S_asw * 10000;

% Display the result
fprintf('\nThe Asset Swap Spread Over Euribor3m is:  %.4f basis points \n', S_asw)

%% Case Study: CDS Bootstrap (Point 2)

% Import data
R = 0.4; % Recovery rate
t0 = dates(1); % Settlement date
% Take the swap dates 1y -> 5y, 7y
datesCDS = [swap1yDate; dates(12:15); dates(17)];
% Spreads in basis points
spreadsCDS = [ 29, 32, 35, 39, 40, 41] / 10000;

%% Complete set of CDS via a spline (Point 2.a)

% Create the complete set of dates
completeDates = [swap1yDate; dates(12:17)];
% Use cubic spline to interpolate the spreads
spreadsCDS = interp1(datesCDS, spreadsCDS, completeDates, 'spline');

% Plot the spreads vs dates
figure
plot(completeDates, spreadsCDS, 'o-')
title('CDS Spreads')
xlabel('Dates')
datetick('x', 'mm/dd/yyyy', 'keepticks')
ylabel('Spreads')

%% Survival probabilities & Intensities (Point 2.b -> 2.d)

% Bootstrap the CDS curve (approx method, neglecting accrual)
[~, P_Approx, int_Approx] = bootstrapCDS(dates, discounts, completeDates, spreadsCDS, 1, R);
% Bootstrap the CDS curve (Exact method, considering accrual)
[~, P_Exact, int_Exact] = bootstrapCDS(dates, discounts, completeDates, spreadsCDS, 2, R);
% Jarrow-Turnbull
[datesCDS, P_JT, int_JT] = bootstrapCDS(dates, discounts, completeDates, spreadsCDS, 3, R);

%(int_Approx + 10^-5) < int_Exact % check and observe that the accrual can
% be negleted, indeed the difference is less than 1bp

% Plot the approx and exact intensities as step functions
% Plot the cumulative mean of the intensities vs JT
PlotIntensities(datesCDS, int_Approx, int_Exact, int_JT)

%% First to Default (Point 3)

% Import data
R_ISP = R;
R_UCG = 0.45;
spreadsCDS_ISP = spreadsCDS;
spreadsCDS_UCG = [34, 39, 45, 46, 47, 47] / 10000; % in bp
rho = 0.2;

% Interpolate the UCG spreads by spline
datesCDS = [swap1yDate; dates(12:15); dates(17)];
spreadsCDS_UCG = interp1(datesCDS, spreadsCDS_UCG, completeDates, 'spline');

% Plot the spreads for ISP and UCG
figure
plot(completeDates, spreadsCDS_ISP, 'o-', 'DisplayName', 'ISP');
hold on;
plot(completeDates, spreadsCDS_UCG, 'o-', 'DisplayName', 'UCG');
title('CDS Spreads');
xlabel('Dates');
datetick('x', 'mm/dd/yyyy', 'keepticks')
ylabel('Spreads');
legend('Location', 'best')

% Take into acocunt the first 4 years
completeDates = completeDates(1:4);
spreadsCDS_ISP = spreadsCDS_ISP(1:4);
spreadsCDS_UCG = spreadsCDS_UCG(1:4);

% Compute the marginal probabilities of survival
% the accrual is neglected, set flag == 2 to have exact calculations
[~, P_ISP, int_ISP] = bootstrapCDS(dates, discounts, completeDates, spreadsCDS_ISP, 1, R_ISP);
[datesCDS, P_UCG, int_UCG] = bootstrapCDS(dates, discounts, completeDates, spreadsCDS_UCG, 1, R_UCG);

% Number of simulations
nSim = 1e4;
% Compute the confidence interval at 95%
confidence_level = 0.95;

% Compute the price of the First to Default (s_FtD)
[priceFtD, lower_bound_FtD, upper_bound_FtD] = priceFirstToDefault(int_ISP, P_ISP, int_UCG, P_UCG, rho, R_ISP,...
    R_UCG, datesCDS, discounts, dates, nSim, confidence_level);

% Display the result
fprintf('\nThe price of the First to Default spread 4 years is :   %.4f bp \n',priceFtD*1e4)
fprintf('Confidence Interval at %d%% confidence level is : [%.4f - %.4f] \n\n', ...
    confidence_level*100, lower_bound_FtD, upper_bound_FtD);

% Different values of the correlation rho
rho = -1:0.20:1;
jj = 0; % counter
h = waitbar(0, 'Calculating prices...');

% Compute the price of the First to Default for each correlation rho
for i = rho
    jj = jj + 1;
    [priceFtD(jj), ~, ~] = priceFirstToDefault(int_ISP, P_ISP, int_UCG, P_UCG, i, ...
        R_ISP, R_UCG, datesCDS, discounts, dates, nSim, confidence_level);
    waitbar(jj / numel(rho), h, sprintf('Progress: %d%%', round(jj / numel(rho) * 100)));
end

close(h);

% Plot the results
figure
plot(rho, priceFtD * 1e4)
title('Price First to Default (Spread)')
xlabel('Correlation (rho)')
ylabel('Spread')

toc