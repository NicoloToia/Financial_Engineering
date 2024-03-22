funtion [VaR_DO, VaR_DM] = graph_report(R)

rng(2);
clc;
clear all;
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

rho = sqrt(R);
% disp('––– Part II Q3: Basel II correlation function –––')
% fprintf('IG correlation: %.3f \n', rho)
% disp(' ')


% Minimum number of Monte Carlo Scenarios
N = 500000;

%% Q1: CreditMetrics' FV (let's start with all bonds rated IG)
FV = FV_risky_bond(IG_cf_schedule_2y, Q, ZC_curve, Recovery); 
E_FV = Q(1,1)*FV(1) + Q(1,2)*FV(2) + Q(1,3)*FV(3);

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



D_count = sum(scen_D)/N;        % Average number of default events in portfolio
d_count = sum(scen_d)/N;        % Average number of "down" events in portfolio
i_count = sum(scen_i)/N;        % Average number of "identical rating" events in portfolio
u_count = sum(scen_u)/N;        % Average number of "up" events in portfolio
% Are MC results in agreement with the input rating transition matrix?
pD_simulated =D_count / N_issuers;
pd_simulated =d_count / N_issuers;
pi_simulated =i_count / N_issuers;
pu_simulated =u_count / N_issuers;

% Default only case (DO)
EL = L_D * pD_simulated + L_d * pd_simulated + L_i * pi_simulated + L_u * pu_simulated;
L_y = quantile(-L_D * scen_D,0.999);
VaR_DO = L_y - EL;

% Default and migration case (DM)
L_y = quantile(-L_D * scen_D - L_d * scen_d - L_i * scen_i,0.999);
VaR_DM = L_y - EL;

