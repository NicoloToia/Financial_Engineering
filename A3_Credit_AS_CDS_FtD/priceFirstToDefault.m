function [priceFtD, lower_bound_FtD, upper_bound_FtD] = priceFirstToDefault(int_ISP, P_ISP,...
    int_UCG, P_UCG, rho, R_ISP, R_UCG, datesCDS, discounts, dates, nSim, confidence_level)
% priceFirstToDefault: price a First to Default (spread_FtD)
%
% INPUT
%   int_ISP         : intensity of the ISP CDS
%   survProbs_ISP   : survival probabilities of the ISP CDS, exp(-int_0&t lambda_ISP(t) dt)
%   int_UCG         : intensity of the UCG CDS
%   survProbs_UCG   : survival probabilities of the UCG CDS, exp(-int_0&t lambda_UCG(t) dt)
%   rho             : correlation between the two CDS
%   R_ISP, R_UCG    : recovery rates for ISP and UCG
%   datesCDS        : dates of the CDS
%   discounts       : discount factors from the Bootstrapping
%   dates           : dates of the Discount Factors
%   nSim            : number of simulations
%   confidence_level: confidence level
%
% OUTPUT
%   priceFtD        : price of the first to default swap
%   lower_bound_FtD : lower bound of the confidence interval
%   upper_bound_FtD : upper bound of the confidence interval

% Compute the discounts at the dates of the CDS
discountsCDS = intExtDF(discounts, dates, datesCDS);
EU_30_360 = 6;
deltas = yearfrac([dates(1); datesCDS(1:end-1)], datesCDS, EU_30_360);

% Simulation via Cholesky decomposition
% A = chol([1 rho; rho 1]);
% y = randn(nSim,2);
% size(A)
% size(y')
% z  = A*y';

% Draw the 2 correlated normal variables (via Matlab function)
z = mvnrnd([0 0], [1 rho; rho 1], nSim);
% Apply the normal cdf to get the correlated uniform variables
u = normcdf(z);
% Create a vector to hold the sum of cashflows of the two legs for each simulation
feeLeg = zeros(nSim, 1);
contingentLeg = zeros(nSim, 1);

completeDates = [dates(1); datesCDS]; % add settlement date

% Invert the two cumulative distribution functions to get the default times
for i=1:nSim
     
    % Invert the survival probabilities to find the default times
    t_ISP = invertProbs(u(i,1), int_ISP, P_ISP, datesCDS);
    t_UCG = invertProbs(u(i,2), int_UCG, P_UCG, datesCDS);

    if isnan(t_ISP) && isnan(t_UCG) % no default happens

        % the default happens after the last date
        % our cashflows are only the full fee leg, not contingent
        % S is omitted, we will compute it later on
        feeLeg(i) = deltas'*discountsCDS;
        % contingentLeg(i) = 0; % not needed, zero by default

%     elseif isnan(t_ISP) % only UCG defaults
% 
%         % only take into consideration the dates before the default time
%         emittedCashflowsIdx = completeDates <= t_UCG;
%         % find the id of the last date before the default time
%         ID_tau = find(emittedCashflowsIdx, 1, 'last');
%         % Find the discount factor at the time of default
%         defaultDF = intExtDF(discounts, dates, t_UCG);
%         % find the fee leg, all payments before the default + the accrual up to default
%         feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) ...
%             + yearfrac(completeDates(ID_tau), t_UCG,EU_30_360)*defaultDF; % accrual term
%         % Find the contingent leg in tau
%         contingentLeg(i) = (1 - R_UCG) * defaultDF;
% 
%     elseif isnan(t_UCG) % only ISP defaults
% 
%         % same as before, only for the ISP rather than the UCG
%         emittedCashflowsIdx = completeDates <= t_ISP;
%         ID_tau = find(emittedCashflowsIdx, 1, 'last');
%         defaultDF = intExtDF(discounts, dates, t_ISP);
%         feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) ...
%             + yearfrac(completeDates(ID_tau), t_ISP, EU_30_360)*defaultDF;
%         % find the contingent leg
%         contingentLeg(i) = (1 - R_ISP) * defaultDF;

    else % both default

        % Find which one happens first
        tau = min(t_ISP, t_UCG);

        emittedCashflowsIdx = completeDates <= tau;
        ID_tau = find(emittedCashflowsIdx, 1, 'last');
        defaultDF = intExtDF(discounts, dates, tau);

        feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) ...
            + yearfrac(completeDates(ID_tau), tau, EU_30_360)*defaultDF;
    
        % Find the contingent leg
        % Find which recovery rate to use
        R = R_ISP*(t_ISP == tau) + R_UCG*(t_UCG == tau);
        contingentLeg(i) = (1 - R)*defaultDF;

        % Check if there is at least one NaN in the contingent leg (security flag)
        if isnan(contingentLeg(i))
            disp('Error: contingent leg is NaN')
            disp(['t_ISP: ' num2str(t_ISP) ' t_UCG: ' num2str(t_UCG) ' tau: ' num2str(tau)])
        end

    end

end

% Write the NPV by taking the expectation of the two legs
% notice: taking the mean of the fee leg written in terms of the ZC bonds is the same 
% as summing the defaultable coupon bonds with the joint CDF.
std_dev_fee = std(feeLeg);  % standard deviation fee leg
feeLeg = mean(feeLeg);      % mean fee leg
% same goes here, taking the mean of the contingent leg is the same as weighing by the probability
% that the default happens in that time tau
std_dev_cont = std(contingentLeg);   % standard devation contingent leg
contingentLeg = mean(contingentLeg); % mean contingent leg

% Finde the price of the First to Default (spread S)
priceFtD = contingentLeg/feeLeg;

% Finde the Confidence Interval
% Given the confidence level compute the quantile
quantile = norminv((1 + confidence_level) / 2, 0, 1);
% Find marginal error fee & contingent leg
error_fee = quantile * (std_dev_fee / sqrt(nSim));
error_cont = quantile * (std_dev_cont / sqrt(nSim));
% Compute endpoints (up/down) fee leg
lower_bound_fee = feeLeg - error_fee;
upper_bound_fee = feeLeg + error_fee;
% Compute endpoints (up/down) contingent leg
lower_bound_cont = contingentLeg - error_cont;
upper_bound_cont = contingentLeg + error_cont;
% Compute endpoints (up/down) First to Default spread
lower_bound_FtD = lower_bound_cont/lower_bound_fee;
upper_bound_FtD = upper_bound_cont/upper_bound_fee;

end