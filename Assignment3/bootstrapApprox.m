function [survProbs, intensities] = bootstrapApprox(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery)
% 

% initialize the output vectors (as column vectors)
survProbs = zeros(length(datesCDS),1);
intensities = zeros(length(datesCDS),1);

% compute the survival probabilities and the intensities
BPV_bar = 0;
sumE = 0;

for i = 1:length(datesCDS)
    % compute the previous probability
    if i == 1
        prevProb = 1;
    else
        prevProb = survProbs(i-1);
    end
    % compute the numerator 
    numerator = (1 - recovery) * (sumE + discountsCDS(i) * prevProb) - BPV_bar * spreadsCDS(i);
    % compute the denominator
    denominator = spreadsCDS(i) * deltas(i) *  discountsCDS(i) + (1 - recovery) * discountsCDS(i);
    % compute the survival probability
    survProbs(i) = numerator / denominator;
    % compute the intensity
    intensities(i) = -log(survProbs(i)/prevProb) / deltasIntensity(i);
    % update the sums
    BPV_bar = BPV_bar + deltas(i) * discounts(i) * survProbs(i);
    sumE = sumE + discounts(i) * (prevProb - survProbs(i));
end

end