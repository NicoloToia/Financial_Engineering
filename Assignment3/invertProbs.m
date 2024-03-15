function tau = invertProbs(u, lambdas, probs, dates)
% Invert the probability of survival to find the time of default
%
% INPUT
%   u        : argument to invert the probability
%   lambdas  : intensities of the survival probability
%   probs    : probabilities of survival
%   dates    : dates of the survival probabilities
% 
% OUTPUT
%   tau      : time to default

% Check if there is default after the last date or before the first date
if u < probs(end)
    tau = NaN;
    return
end

% Find the index of the following date
nextIdx = find(u>probs, 1, 'first');

% Find the exact tau by inverting the exponential and going back
fraction = -log(probs(nextIdx)/u)/lambdas(nextIdx);
% Invert the year frac
% (Order is very important!!!, otherwise we get dates after the next date instead of previous dates
% as we have assumed)
ACT_365 = 3;
fun = @(x) yearfrac(x, dates(nextIdx), ACT_365) - fraction;
tau = fzero(fun, dates(nextIdx));
% Round up
tau = ceil(tau);

% Check that tau is less than the next date
%disp(['tau = ', datestr(tau, 'dd/mm/yyyy'), ' and next date = ', datestr(dates(nextIdx), 'dd/mm/yyyy')]);

end