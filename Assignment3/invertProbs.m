function tau = invertProbs(u, lambdas, probs, dates)
% invert the probability of survival to find the time of default
%
% Inputs
% u - argument to invert the probability
% lambdas - intensities of the survival probability
% probs - probabilities of survival
% dates - dates of the survival probabilities

% find the probability that is just less than u
nextIdx = find(probs < u, 1, 'last');

% if the index is more than the length of the probabilities
% return none
if isempty(nextIdx)
    tau = NaN;
    return
end

% find the exact tau by inverting the exponential and going back
fraction = -log(probs(nextIdx)/u)/lambdas(nextIdx);
% invert the year frac
% (order is very important!!!, otherwise we get dates after the next date instead of previous dates
% as we have assumed)
fun = @(x) yearfrac(x, dates(nextIdx)) - fraction;
tau = fzero(fun, dates(nextIdx));

% convert to a date
tau = datenum(tau);

end