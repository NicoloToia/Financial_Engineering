function optionPrice = BermudanOptionCRR_new(F0, K, B, T, sigma, N, flag)
% Bermudan option price with CRR method
%
% INPUT
% F0:    forward price
% K:     strike
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% flag:  1 call, -1 put
%
% OUTPUT
% optionPrice : Price of the Bermudan option

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

% initialize the tree leaves (N+1) with the payoff
Ftt = F0 * u.^(N:-2:-N);
leavesCRR = max(Ftt - K, 0);
S0 = 1;

% reduce the tree to the root
for i = 1:N
    if i == (round(N/3) + 1) || i == round((2/3)*N) + 1
        % vectorized version
        ST_t = max(S0 *u.^(N-i:-2:-N+i) - K,0);
        leavesCRR = max(B_dt *(q * leavesCRR(1:end-1) + (1-q) * leavesCRR(2:end)),ST_t);
    else
        leavesCRR = B_dt *(q * leavesCRR(1:end-1) + (1-q) * leavesCRR(2:end));
    end
end

CallPrice = leavesCRR;

% exploit put-call parity
if flag == 1
    optionPrice = CallPrice;
else
    optionPrice = CallPrice - B * (F0 - K);
end

end % function EuropeanOptionCRR