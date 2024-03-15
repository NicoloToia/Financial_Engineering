function [survProbs, intensities] = bootstrapExact(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery)
% This function computes the survival probabilities and the intensities
% considering the accrual term in the CDS contract.
%
% INPUTS
%   discountsCDS    : a vector of discounts at the CDS payment dates
%   spreadsCDS      : a vector of CDS spreads
%   deltas          : a vector of time fractions between the CDS payment dates (EU_30_360 convention)
%   deltasIntensity : same as deltas but with the ACT/365 convention
%   recovery        : the recovery rate
% OUTPUTS
%   survProbs       : a vector of survival probabilities
%   intensities     : a vector of intensities(lambda)

% Initialize the output vectors (as column vectors)
survProbs = zeros(length(spreadsCDS),1);
intensities = zeros(length(spreadsCDS),1);

% Compute the survival probabilities and the intensities
BPV_bar = 0; % sum of B(t0;t_i) * P(t_i) up to n-1
sumE = 0; % sum of the e(t0;t_{i-1},t_i) up to n-1
sumAccrual = 0; % sum of delta_i * e(t0;t_{i-1},t_i) up to n-1

for n = 1:length(spreadsCDS)
    if n == 1
        prevProb = 1;
    else
        prevProb = survProbs(n-1);
    end
    % Compute the numerator 
    numerator = (1 - recovery) * (sumE + discountsCDS(n) * prevProb) - ...
        spreadsCDS(n) * ( BPV_bar + sumAccrual + deltas(n) /2 * discountsCDS(n) * prevProb);
    % Compute the denominator
    denominator = spreadsCDS(n) * deltas(n) / 2 *  discountsCDS(n) + (1 - recovery) * discountsCDS(n);
    % Compute the survival probability
    survProbs(n) = numerator / denominator;
    % compute the intensity
    intensities(n) = -log(survProbs(n)/prevProb) / deltasIntensity(n);
    % Update the sums
    BPV_bar = BPV_bar + deltas(n) * discountsCDS(n) * survProbs(n);
    sumE = sumE + discountsCDS(n) * (prevProb - survProbs(n));
    sumAccrual = sumAccrual + deltas(n) / 2 * (prevProb - survProbs(n));
end

end