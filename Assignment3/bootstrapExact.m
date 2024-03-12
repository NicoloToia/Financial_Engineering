function [survProbs, intensities] = bootstrapExact(discountsCDS, spreadsCDS, deltas, deltasIntensity, recovery)

% initialize the output vectors (as column vectors)
survProbs = zeros(length(spreadsCDS),1);
intensities = zeros(length(spreadsCDS),1);

% compute the survival probabilities and the intensities
BPV_bar = 0; % sum of B(t0;t_i) * P(t_i) up to n-1
sumE = 0; % sum of the e(t0;t_{i-1},t_i) up to n-1
sumAccrual = 0; % sum of delta_i * e(t0;t_{i-1},t_i) up to n-1

for n = 1:length(spreadsCDS)
    if n == 1
        prevProb = 1;
    else
        prevProb = survProbs(n-1);
    end
    % compute the numerator 
    numerator = (1 - recovery) * (sumE + discountsCDS(n) * prevProb) - ...
        spreadsCDS(n) * ( BPV_bar + sumAccrual + deltas(n) /2 * discountsCDS(n) * prevProb);
    % compute the denominator
    denominator = spreadsCDS(n) * deltas(n) / 2 *  discountsCDS(n) + (1 - recovery) * discountsCDS(n);
    % compute the survival probability
    survProbs(n) = numerator / denominator;
    % compute the intensity
    intensities(n) = -log(survProbs(n)/prevProb) / deltasIntensity(n);
    % update the sums
    BPV_bar = BPV_bar + deltas(n) * discountsCDS(n) * survProbs(n);
    sumE = sumE + discountsCDS(n) * (prevProb - survProbs(n));
    sumAccrual = sumAccrual + deltas(n) / 2 * (prevProb - survProbs(n));
end

end