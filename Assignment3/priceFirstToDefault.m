function [priceFtD, lower_bound_FtD, upper_bound_FtD] = priceFirstToDefault(int_ISP, P_ISP,...
    int_UCG, P_UCG, rho, R_ISP, R_UCG, datesCDS, discounts, dates, nSim, confidence_level)
% priceFirstToDefault: price a first to default swap
%
% INPUT:
% int_ISP: intensity of the ISP CDS
% survProbs_ISP: survival probabilities of the ISP CDS, exp(-int_0&t lambda_ISP(t) dt)
% int_UCG: intensity of the UCG CDS
% survProbs_UCG: survival probabilities of the UCG CDS, exp(-int_0&t lambda_UCG(t) dt)
% rho: correlation between the two CDS
% R_ISP, R_UCG: recovery rates for ISP and UCG
% datesCDS: dates of the CDS
% discounts: discount factors from the Bootstrapping
% dates: dates of the Discount Factors
%
% OUTPUT:
% priceFtD: price of the first to default swap

% compute the discounts at the dates of the CDS
discountsCDS = intExtDF(discounts, dates, datesCDS);
EU_30_360 = 6;
deltas = yearfrac([dates(1); datesCDS(1:end-1)], datesCDS, EU_30_360);

% Simulation via Cholesky decomposition
% A = chol([1 rho; rho 1]);
% y = randn(nSim,2);
% size(A)
% size(y')
% z  = A*y';

% draw the 2 correlated normal variables (via Matlab function)
z = mvnrnd([0 0], [1 rho; rho 1], nSim);
% apply the normal cdf to get the correlated uniform variables
u = normcdf(z);
% create a vector to hold the sum of cashflows of the two legs for each simulation
feeLeg = zeros(nSim, 1);
contingentLeg = zeros(nSim, 1);

completeDates = [dates(1); datesCDS];

% invert the two cumulative distribution functions to get the default times
nDefaults = 0;
for i=1:nSim
     
    % invert the survival probabilities to find the default times
    t_ISP = invertProbs(u(i,1), int_ISP, P_ISP, datesCDS);
    t_UCG = invertProbs(u(i,2), int_UCG, P_UCG, datesCDS);

    % count how many defaults we have
    nDefaults = nDefaults + sum(~isnan([t_ISP t_UCG]));

    if isnan(t_ISP) && isnan(t_UCG) % no default happens

        % the default happens after the last date
        % our cashflows are only the full fee leg, not contingent
        % S is omitted, we will compute it later on
        feeLeg(i) = deltas'*discountsCDS;
        contingentLeg(i) = 0;                                               % ADD CONTINGENT ZERO INUTILE è GIA ZERO MA PER CAPISSE

    elseif isnan(t_ISP) % only UCG defaults

%         % find the cashflows (as the payments done before the default) of the fee leg
%         emittedCashflowsIdx = datesCDS <= t_UCG;
%         %emittedCashflowsIdx(find(emittedCashflowsIdx == 0, 1, 'first')) = 1;
%         feeLeg(i) = deltas(emittedCashflowsIdx)'*discountsCDS(emittedCashflowsIdx);

        % find the cashflows (as the payments done before the default) of the fee leg
        emittedCashflowsIdx = completeDates <= t_UCG;
        ID_tau = find(emittedCashflowsIdx, 1, 'last');
        defaultDF = intExtDF(discounts, dates, t_UCG);
        feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) + yearfrac(completeDates(ID_tau), t_UCG,EU_30_360)*defaultDF;
        % find the contingent leg
        contingentLeg(i) = (1 - R_UCG)*defaultDF;

%         % find the contingent leg
%         defaultDF = discountsCDS(find(datesCDS >= t_UCG, 1, 'first'));
%         contingentLeg(i) = (1 - R_UCG)*defaultDF;

    elseif isnan(t_UCG) % only ISP defaults

%         % find the cashflows before
%         emittedCashflowsIdx = datesCDS <= t_ISP;
%         %emittedCashflowsIdx(find(emittedCashflowsIdx == 0, 1, 'first')) = 1;
%         feeLeg(i) = deltas(emittedCashflowsIdx)'*discountsCDS(emittedCashflowsIdx);
% 
%         defaultDF = discountsCDS(find(datesCDS >= t_ISP, 1, 'first'));
%         contingentLeg(i) = (1 - R_ISP)*defaultDF;

        emittedCashflowsIdx = completeDates <= t_ISP;
        ID_tau = find(emittedCashflowsIdx, 1, 'last');
        defaultDF = intExtDF(discounts, dates, t_ISP);
        feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) + yearfrac(completeDates(ID_tau), t_ISP, EU_30_360)*defaultDF;
        % find the contingent leg
        contingentLeg(i) = (1 - R_ISP)*defaultDF;

    else % both default

        % find which one happens first
        tau = min(t_ISP, t_UCG);

        % find the cashflows (as the payments done before the default)
        emittedCashflowsIdx = completeDates <= tau;
        ID_tau = find(emittedCashflowsIdx, 1, 'last');

        defaultDF = intExtDF(discounts, dates, tau);
        %emittedCashflowsIdx(find(emittedCashflowsIdx == 0, 1, 'first')) = 1;
        feeLeg(i) = deltas(emittedCashflowsIdx(2:end))'*discountsCDS(emittedCashflowsIdx(2:end)) + yearfrac(completeDates(ID_tau), tau, EU_30_360)*defaultDF;
    
        % find the contingent leg
        % find which recovery rate to use
        R = R_ISP*(t_ISP == tau) + R_UCG*(t_UCG == tau);
        contingentLeg(i) = (1 - R)*defaultDF;

        % check if there is at least one NaN in the contingent leg
        if isnan(contingentLeg(i))
            disp('Error: contingent leg is NaN')
            disp(['t_ISP: ' num2str(t_ISP) ' t_UCG: ' num2str(t_UCG) ' tau: ' num2str(tau)])
        end

    end

end

%     priceFtD = contingentLeg./feeLeg;               % PROVA DA LEVARE
%     priceFtD(isinf(priceFtD)) = 0;   
%     
%     priceFtD = mean(priceFtD);               % ERRORE

% write the NPV by taking the expectation of the two legs
% notice: taking the mean of the fee leg written in terms of the ZC bonds is the same 
% as summing the defaultable coupon bonds with the joint CDF.

std_dev_fee = std(feeLeg);
feeLeg = mean(feeLeg);
% same goes here, taking the mean of the contingent leg is the same as weighing by the probability
% that the default happens in that time tau
std_dev_cont = std(contingentLeg);
contingentLeg = mean(contingentLeg);

priceFtD = contingentLeg/feeLeg;

quantile = norminv((1 + confidence_level) / 2, 0, 1);

margin_error_fee = quantile * (std_dev_fee / sqrt(nSim));
margin_error_cont = quantile * (std_dev_cont / sqrt(nSim));

lower_bound_fee = feeLeg - margin_error_fee;
upper_bound_fee = feeLeg + margin_error_fee;

lower_bound_cont = contingentLeg - margin_error_cont;
upper_bound_cont = contingentLeg + margin_error_cont;

lower_bound_FtD = lower_bound_cont/lower_bound_fee;
upper_bound_FtD = upper_bound_cont/upper_bound_fee;

end