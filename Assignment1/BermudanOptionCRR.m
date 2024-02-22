function OptionPrice = BermudanOptionCRR(F0, K, B, T, sigma, N, flag)
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
%q = (1-1/u) / (u-1/u);
q = 1 / (u + 1);

r = -log(B) / T;
B_dt = exp(-r*dt);

% initialize the tree leaves (N+1) with the payoff
Ftt = F0 * u.^(N:-2:-N);
leavesCRR = max(Ftt - K, 0);

% reduce the tree to the root
for i = 1:N
    % vectorized version
    leavesCRR = max(B_dt *(q * leavesCRR(1:end-1) + (1-q) * leavesCRR(2:end)),max());
end
CallPrice = leavesCRR;

% exploit put-call parity
if flag == 1
    OptionPrice = CallPrice;
else
    OptionPrice = CallPrice - B * (F0 - K);
end

end % function EuropeanOptionCRR