function payoffplot = plotpayoff(strikes, k)
% Plot the payoff graph for a digital call option
% INPUT:
% strikes: vector of strike prices
% k: ATM strike price
% OUTPUT:
% payoffplot: plot for payoff graph for a digital call option and a BUll spread


% Parameters
S = 2600:0.1:strikes(end); % Underlying asset's price range
K = k; % Strike price
payoff_digital_call = (S > K); % Payoff for digital call option

% Find the index where S exceeds K
idx = find(S > K, 1);

% Plot
figure;
plot(S(1:idx-1), zeros(size(S(1:idx-1))), 'r', 'LineWidth', 2); % Payoff is 0 before K
hold on;
plot(S(idx:end), ones(size(S(idx:end))), 'r', 'LineWidth', 2); % Payoff is 1 after K
xlabel('Underlying Asset Price');
ylabel('Payoff');
title('Payoff Graph for Digital Call Option');
grid on;


% Add filled and empty circles
plot(S(idx), 1, 'ro', 'MarkerFaceColor', 'r'); % Filled circle at jump point
plot(S(idx-1), 0, 'ro', 'MarkerFaceColor', 'w'); % Empty circle at jump point

ylim([0, 1.2]);

% Insert K in the x-axis
xticks([k]);
xticklabels({['K']});

% Parameters
delta_k = 200; % Strike price increment
K_lower = k; % Lower strike price
K_upper = k + delta_k; % Upper strike price
payoff_bull_spread = max(min(S - K_lower, 1), 0) - max(min(S - K_upper, 1), 0); % Payoff for bull call spread

% Payoff for long call option
payoff_long_call = max(S - K_lower, 0);
% Payoff for short call option
payoff_short_call = max(S - K_upper, 0);

% Payoff for bull call spread
payoff_bull_spread = payoff_long_call - payoff_short_call;
payoff_bull_spread = payoff_bull_spread/delta_k;

% Plot
figure;
plot(S, payoff_bull_spread, 'b', 'LineWidth', 2); % Payoff for bull call spread
xlabel('Underlying Asset Price');
ylabel('Payoff');
title('Payoff Graph for Bull Call Spread');
grid on;

% Set y-axis limits
ylim([0, 1.2]);

% Add x-axis ticks for K and K + delta_K
xticks([K_lower, K_upper]);
xticklabels({['K'], ['K + \Delta K']});







payoffplot = 1;
