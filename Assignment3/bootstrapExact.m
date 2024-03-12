function [survProbs, intensities] = bootstrapExact(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery)

% initialize the output vectors (as column vectors)
survProbs = zeros(length(spreadsCDS),1);
intensities = zeros(length(spreadsCDS),1);

% compute the survival probabilities and the intensities
BPV_bar = 0; % sum of B(t0;t_i) * P(t_i) up to n-1
sumE = 0; % sum of the e(t0;t_{i-1},t_i) up to n-1
sumAccrual = 0; 

for i = 1:length(spreadsCDS)
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