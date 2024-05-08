function [price_digital] = DigitalPrice(Notional , T , F_0 , dividend, discount_1y , k , strikes , surface , flag)

% Digital_Price: Computes the price of a digital option via different
% approaches
%
% INPUT: 
% Notional      : Notional of the option
% T             : Time to maturity
% F_0           : Forward price
% discount_1y   : Discount factor at 1 year
% k             : Strike price
% strikes       : Vector of strikes
% surface       : Vector of volatilities (smile)
% flag          : Flag to choose the pricing method
%                   - flag = 1 --> Black formula 
%                   - flag = 2 --> Implied Volatility Approach                      
%                   - flag = 3 --> Monte Carlo (Black dynamics)

% OUTPUT:
% price_digital : Price of the digital option

% find the corresponding volatility from the smile
sigma_digital = interp1(strikes, surface, k, 'spline');

% value of the digital option at maturity if s > k
payment = Notional;

% compute d_1 and d_2 Black's formula
d_1 = (log(F_0 / k) + (0.5 * sigma_digital^2) * T) / (sigma_digital * sqrt(T));
d_2 = d_1 - sigma_digital * sqrt(T);

if flag==1 % Black Formula

    % The closed formula under Black
    price_digital = payment * discount_1y * normcdf(d_2);

elseif flag==2 % Implied Volatility Approach

    % compute the skew for the entire surface
    skew = zeros(size(strikes));
    skew(1) = (surface(2) - surface(1)) / (strikes(2) - strikes(1));
    
    for i = 2:length(strikes) - 1
        skew(i) = (surface(i + 1) - surface(i - 1)) / (strikes(i + 1) - strikes(i - 1));
    end

    % plot the skew vs the strikes
    % figure;
    % plot(strikes, skew);
    % xlabel('Strikes'); 
    % ylabel('Skew');
    % title('Skew vs Strikes');
    
    % find the skew in the given strike
    m = interp1(strikes, skew, k, 'spline');

    % Compute the vega under black model
    d_1 = (log(F_0 ./ strikes) + (0.5 * surface.^2) * T) ./ (surface .* sqrt(T));
    vega = F_0 .* discount_1y .* normpdf(d_1) * sqrt(T);

    % plot the vega vs the strikes
    % figure;
    % plot(strikes, vega);
    % xlabel('Strikes');
    % ylabel('Vega');
    % title('Vega vs Strikes');

    % find the vega in the given strike
    vega_k = interp1(strikes, vega, F_0, 'spline');

    % compute the price of the digital option under Black
    price_digital_black = payment * discount_1y * normcdf(d_2);

    % Now compute the digital price with the implied volatility correction
    price_digital = price_digital_black - vega_k * m * payment;

elseif flag==3 % MC Black dynamics

    % now implement monte carlo simulation to check whether the price under Black is correct
    N = 1e7;
    Z = randn(N, 1);
    F_t = F_0 * exp(-0.5 * sigma_digital^2 * T + sigma_digital * sqrt(T) * Z);
    % payoff digital option pays 0.05 
    payoff = payment *  (F_t > k);
    price_digital = mean(payoff) * discount_1y;

else
    disp('Flag not recognized');
end

end