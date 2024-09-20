function optionPrice = EuropeanOptionCRR(F0, K, B, T, sigma, N, flag)
%European option price with CRR method
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N :    number of time steps (knots for CRR tree)
% flag:  1 call, -1 put
%
% OUTPUT
% optionPrice : Price of the option

% Compute the following parameters:
% dt is the interval of time between two knots
dt = T/N;
% u is the factor by which the stock price increases
u = exp(sigma * sqrt(dt));
% q is the probability under the risk neutral measure
q = 1 / (u + 1);

% Compute the discount factor for each step
r = -log(B) / T;
B_dt = exp(-r*dt);

% Initialize the tree leaves (N+1) with the payoff
Ftt = F0 * u.^(N:-2:-N);
leavesCRR = max(Ftt - K, 0);

% Reduce the tree to the root to find the option price
for i = 1:N
    % vectorized version
    leavesCRR = B_dt * (q * leavesCRR(1:end-1) + (1-q) * leavesCRR(2:end));
end

CallPrice = leavesCRR;

% Exploit put-call parity: C0 - P0 = B*(F0 - K)
if flag == 1
    optionPrice = CallPrice;                    % Call price
else
    optionPrice = CallPrice - B * (F0 - K);     % Put price
end

end % function EuropeanOptionCRR