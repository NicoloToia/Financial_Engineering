function [datesCDS, survProbs, intensities] =  bootstrapCDS_v2(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)

% INPUT:   
%   datesDF: dates of the discount factors
%   discounts: discount factors
%   datesCDS: dates of the CDS
%   spreadsCDS: spreads of the CDS
%   flag: 1 (approx), 2 (exact), 3 (Jarrow-Turnbull)
%   recovery: recovery rate

% OUTPUT: 
%   datesCDS: dates of the CDS
%   survProbs: survival probabilities
%   intensities: intensitiesz


% compute the discount factors at the CDS dates
discountsCDS = interp1(datesDF,discounts,datesCDS);
% compute the year fractions
EU_30_360 = 6;
ACT_365 = 3;
deltas = yearfrac([datesDF(1); datesCDS(1:end-1)], datesCDS, EU_30_360);
deltasIntensity = yearfrac([datesDF(1); datesCDS(1:end-1)], datesCDS, ACT_365);

% switch the flag
switch (flag)
    case 1 % approximate (neglect accrual)
        [survProbs, intensities] = bootstrapApprox(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery);
    case 2 % exact (consider accrual)
        [survProbs, intensities] = bootstrapExact(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery);
    case 3 % Jarrow-Turnbull
        [survProbs, intensities] = bootstrapJT(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery);
    otherwise
        % throw an error
        error('Flag not supported');
end


end