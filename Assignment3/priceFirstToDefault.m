function priceFtD = priceFirstToDefault(int_ISP, P_ISP, int_UCG, P_UCG, rho, R_ISP, R_UCG, datesCDS, discounts, dates)
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

% number of simulations
nSim = 1000;

% draw the 2 correlated normal variables
z = mvnrnd([0 0], [1 rho; rho 1], nSim);
% apply the normal cdf to get the correlated uniform variables
u = normcdf(z);
% create a vector to hold the sum of cashflows of the two legs for each simulation
feeLeg = zeros(nSim, 1);
contingentLeg = zeros(nSim, 1);

% invert the two cumulative distribution functions to get the default times
for i=1:nSim
     
    % invert the survival probabilities to find the default times
    t_ISP = invertProbs(u(i,1), int_ISP, P_ISP, datesCDS);
    t_UCG = invertProbs(u(i,2), int_UCG, P_UCG, datesCDS);

    % check that at least one of them is not NaN
    if isnan(t_ISP) && isnan(t_UCG)

        % the default happens after the last date
        % our cashflows are only the full fee leg, not contingent
        % S is omitted, we will compute it later on
        feeLeg(i) = deltas'*discountsCDS;

    elseif isnan(t_ISP) % only UCG defaults

        % find the cashflows (as the payments done before the default) of the fee leg
        emittedCashflowsIdx = datesCDS < t_UCG;
        feeLeg(i) = deltas(emittedCashflowsIdx)'*discountsCDS(emittedCashflowsIdx);
        % find the contingent leg
        defaultDF = intExtDF(discounts, dates, t_UCG);
        contingentLeg(i) = (1 - R_UCG)*defaultDF;

    elseif isnan(t_UCG) % only ISP defaults

        emittedCashflowsIdx = datesCDS < t_ISP;
        feeLeg(i) = deltas(emittedCashflowsIdx)'*discountsCDS(emittedCashflowsIdx);

        defaultDF = intExtDF(discounts, dates, t_ISP);
        contingentLeg(i) = (1 - R_ISP)*defaultDF;

    else % both default

        % find which one happens first
        tau = min(t_ISP, t_UCG);

        % find the cashflows (as the payments done before the default)
        emittedCashflowsIdx = datesCDS < tau;
        feeLeg(i) = deltas(emittedCashflowsIdx)'*discountsCDS(emittedCashflowsIdx);

        % find the contingent leg
        defaultDF = intExtDF(discounts, dates, tau);
        % find which recovery rate to use
        R = R_ISP*(t_ISP == tau) + R_UCG*(t_UCG == tau);
        contingentLeg(i) = (1 - R)*defaultDF;

    end

end

% write the NPV by taking the expectation of the two legs
feeLeg = mean(feeLeg);
contingentLeg = mean(contingentLeg);
NPV = @(S) S * feeLeg - contingentLeg;

% neutralize the NPV to find the spread
priceFtD = fzero(NPV, 0.5);

end