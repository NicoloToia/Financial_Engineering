function [price_digital] = Digital_Price(Notional , T , F_0 , discount_1y , sigma_digital , k , strikes , surface , flag);
% Digital_Price: Computes the price of a digital option
%
% INPUT: 
% Notional: Notional of the option
% T: Time to maturity
% F_0: Forward price
% discount_1y: Discount factor at 1 year
% sigma_digital: Volatility of the digital option
% k: Strike price
% strikes: Vector of strikes
% surface: Vector of volatilities
% flag: Flag to choose the pricing method

% OUTPUT:
% price_digital: Price of the digital option

% value of the digital option at maturity if s > k
payment = 0.05 * Notional;
d_1 = (log(F_0 / k) + (0.5 * sigma_digital^2) * T) / (sigma_digital * sqrt(T));
d_2 = d_1 - sigma_digital * sqrt(T);

if flag==1
    price_digital = payment * discount_1y * normcdf(d_2);
end
if flag==2
    % Now use a volatility approach
    % Find k_1 and k_2 that contain k
    k_1 = strikes(find(strikes < k, 1, 'last'));
    k_2 = strikes(find(strikes > k, 1, 'first'));

    % Find the corresponding volatilities
    sigma_1 = interp1(strikes, surface, k_1, 'spline');
    sigma_2 = interp1(strikes, surface, k_2, 'spline');

    % compute the skew in that point
    m = (sigma_2 - sigma_1) / (k_2 - k_1);

    % Compute the vega under black model
    vega = F_0 * discount_1y * normpdf(d_1) * sqrt(T) * 0.01;

    price_digital_black = payment * discount_1y * normcdf(d_2);

    % Now compute the digital price
    price_digital = price_digital_black - vega * m * payment;
end
if flag==3
    % now implement monte carlo simulation
    N = 1e7;
    Z = randn(N, 1);
    F_t = F_0 * exp(-0.5 * sigma_digital^2 * T + sigma_digital * sqrt(T) * Z);
    % payoff digital option pays 0.05 
    payoff = payment *  (F_t > k);
    price_digital = mean(payoff) * discount_1y;
end
