function [price_digital] = Digital_Price(Notional , T , F_0 , dividend, discount_1y , sigma_digital , k , strikes , surface , flag);
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

    % compute the skew for the entire surface
    skew = zeros(size(strikes));
    skew(1) = (surface(2) - surface(1)) / (strikes(2) - strikes(1));
    
    for i = 2:length(strikes) - 1
        skew(i) = (surface(i + 1) - surface(i - 1)) / (strikes(i + 1) - strikes(i - 1));
    end

    % plot the skew vs the strikes
    figure;
    plot(strikes, skew);
    xlabel('Strikes'); 
    ylabel('Skew');
    title('Skew vs Strikes');
    
    % find the skew in that point
    m = interp1(strikes, skew, k, 'spline');

    % Compute the vega under black model
    F0_values = strikes/discount_1y * exp(-dividend * T);
    d_1 = (log(F0_values ./ k) + (0.5 * surface.^2) * T) ./ (surface .* sqrt(T));
    vega = F0_values .* discount_1y .* normpdf(d_1) * sqrt(T).* 0.01;

    % plot the vega vs the strikes
    figure;
    plot(strikes, vega);
    xlabel('Strikes');
    ylabel('Vega');
    title('Vega vs Strikes');

    vega_k = interp1(strikes, vega, F_0, 'spline');

    price_digital_black = payment * discount_1y * normcdf(d_2);

    % Now compute the digital price
    price_digital = price_digital_black - vega_k * m * payment;
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
