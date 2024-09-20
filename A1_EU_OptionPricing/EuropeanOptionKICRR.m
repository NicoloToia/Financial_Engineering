function optionPrice=EuropeanOptionKICRR(F0,K,KI,B,T,sigma,N)
% EuropeanOptionKICRR computes the price of a European option with a
% knock-in barrier using the CRR method
%
% INPUT
% F0:    forward price
% K:     strike
% KI:    barrier
% B:     discount factor
% T:     time-to-maturity
% sigma: volatility
% N:     number time steps (knots for CRR tree)
%
% OUTPUT
% optionPrice : Price of the option with a knock-in barrier

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

% Initialize the tree leaves with the payoff of the call option with
% European barrier KI
Ftt = F0 * u.^(N:-2:-N);
leavesCRR = max((Ftt - K), 0) .* (Ftt>KI);

% Reduce the tree to the root to find the option price
for i = 1:N
    % vectorized version
    leavesCRR = B_dt * (q*leavesCRR(1:end-1) + (1-q)*leavesCRR(2:end));
end

% Return the option price
optionPrice = leavesCRR;

end