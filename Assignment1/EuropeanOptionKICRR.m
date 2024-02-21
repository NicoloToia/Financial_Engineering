function optionPrice = EuropeanOptionKICRR(F0,K, KI,B,T,sigma,N)
%INPUT
% F0:    forward price
% B:     discount factor
% K:     strike
% KI:    barrier
% T:     time-to-maturity
% sigma: volatility
% N:     either number of time steps (knots for CRR tree)

% compute the parameters
dt = T/N;
u = exp(sigma * sqrt(dt));
d = 1/u;
q = (1-d) / (u-d);

% compute the value of the option at the leaves and the probability of each leaf
% for N steps, there are N+1 leaves
CRR_leaves = zeros(N+1, 1);
CRR_prob = zeros(N+1, 1);

% m is the number of up moves
for m = 0:N
    % compute the value of the forward at the leaf
    Ftt = F0 * u^(m) * d^(N-m);
    % compute the option value
    CRR_leaves(m+1) = B*max((Ftt - K), 0)*(Ftt>KI);
    % compute the probability of this leaf
    % the probability of a leaf is (N choose m) * q^m * (1-q)^(N-m)
    CRR_prob(m+1) = nchoosek(N, m) * q^m * (1-q)^(N-m);
end

% compute the option price
% Price is the mean of the leaves (weighted by their probability) discounted
C0 = sum(CRR_leaves .* CRR_prob);

if flag == 1
    optionPrice = C0;
else % put option
    % leverage the put-call parity
    % C0 - P0 = B*(F0 - K)
    optionPrice = C0 - B*(F0 - K);
end
end
