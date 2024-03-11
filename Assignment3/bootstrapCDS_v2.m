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

% initialize the output vectors (as column vectors)
survProbs = zeros(length(datesCDS),1);
intensities = zeros(length(datesCDS),1);

% compute the discount factors at the CDS dates
discountsCDS = interp1(datesDF,discounts,datesCDS);
% compute the year fractions
EU_30_360 = 6;
deltas = yearfrac([datesDF(1); datesCDS(1:end-1)], datesCDS, EU_30_360);

% compute the survival probabilities and the intensities
BPV_bar = 0;
sumE = 0;

for i = 1:length(datesCDS)
    if i == 1
        prevProb = 1;
    else
        prevProb = survProbs(i-1);
    end
    % compute the numerator 
    numerator = (1 - recovery) * (sumE + discountsCDS(i) * prevProb) - BPV_bar * spreadsCDS(i);
    % compute the denominator
    denominator = spreadsCDS(i) * deltas(i) *  discountsCDS(i) + (1 - recovery) * discountsCDS(i);
    survProbs(i) = numerator / denominator;
    % update the sums
    BPV_bar = BPV_bar + deltas(i) * discounts(i) * survProbs(i);
    sumE = sumE + discounts(i) * (prevProb - survProbs(i));
end

end