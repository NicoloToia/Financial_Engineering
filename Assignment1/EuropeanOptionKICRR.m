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
q = (1-1/u) / (u-1/u);

r = -log(B) / T;
B_dt = exp(-r*dt);

% initialize the tree leaves with the payoff
Ftt = F0 * u.^(N:-2:-N);
leavesCRR = max((Ftt - K), 0) .* (Ftt>KI);

% go back to the root
for i = 1:N
    leavesCRR = B_dt * (q*leavesCRR(1:end-1) + (1-q)*leavesCRR(2:end));
end

% return the option price
optionPrice = leavesCRR;

end
