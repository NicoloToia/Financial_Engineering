function [survProbs, intensities] = bootstrapApprox(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery)
% This function computes the survival probabilities and the intensities
% ignoring the accrual term in the CDS contract.
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
BPV_bar = 0; % sum of B(t_0, t_i) * P(t_0, t_i) up to n-1
sumE = 0; % sum of E(t_0, t_{i-1}, t_i) up to n-1

% For each CDS payment date
for n = 1:length(spreadsCDS)
    % Compute the previous probability (by definition P(t_0, t_0) = 1)
    if n == 1
        prevProb = 1;
    else
        prevProb = survProbs(n-1);
    end
    % Compute the numerator 
    % (1 - R) * (sumE + B(t_0, t_n) * P(t_0, t_n)) - BPV_var * S(t_0, t_n)
    numerator = (1 - recovery) * (sumE + discountsCDS(n) * prevProb) - BPV_bar * spreadsCDS(n);
    % Compute the denominator
    % S(t_0, t_n) * B(t_0, t_n) * delta(t_{i-1}, t_i) + (1 - R) * B(t_0, t_n)
    denominator = spreadsCDS(n) * deltas(n) *  discountsCDS(n) + (1 - recovery) * discountsCDS(n);
    % Compute the survival probability
    survProbs(n) = numerator / denominator;
    % Compute the intensity
    intensities(n) = -log(survProbs(n)/prevProb) / deltasIntensity(n);
    % Update the sums
    BPV_bar = BPV_bar + deltas(n) * discountsCDS(n) * survProbs(n);
    sumE = sumE + discountsCDS(n) * (prevProb - survProbs(n));
end

end