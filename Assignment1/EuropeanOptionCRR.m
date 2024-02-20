function OptionPrice = EuropeanOptionCRR(F0, K, B, T, sigma, N, flag)
%European option price with CRR method
%
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% T:     time-to-maturity
% sigma: volatility
% flag:  1 call, -1 put

% compute the parameters
dt = T/N;
u = exp(sigma * sqrt(dt));
d = 1/u;
q = (1-d) / (u-d);

% compute the value of the forward at the leaves
leaves = zeros(N+1, 1);

for i = 0:N
    leaves(i+1) = F0 * u^(i) * d^(N-i);
end

% compute the mean of the payoff at the leaves
% m represents the number of up movements
for m = 0:N
    % compute the option value
    optionValue = max(flag * (leaves(m+1) - K), 0);
    % compute the probability of this leaf
    % the probability of a leaf is (N choose m) * q^m * (1-q)^(N-m)
    prob = nchoosek(N, m) * q^m * (1-q)^(N-m);
    % store the value
    leaves(m+1) = optionValue * prob;
end

% compute the option price
% Price is the mean of the leaves discounted
OptionPrice = B * sum(leaves);

end