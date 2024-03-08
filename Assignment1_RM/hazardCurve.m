function [h_curve] = hazardCurve(zeroCurve, R, dirtyPrice_1y, dirtyPrice_2y, ...
    couponSchedule_1y, couponSchedule_2y)
% Bootstrap the hazard 
% 
% INPUT
% zeroCurve : ZC bond data [yearfrac ; rates]
% R : recovery rate
% dirtyPrice_1y : corporate market price 1year bond
% dirtyPrice_2y : corporate market price 2year bond
% couponSchedule_1y : cash flows 1year bond [yearfrac; cash flow]
% couponSchedule_2y : cash flows 2year bond [yearfrac; cash flow]
%
% OUTPUT
% h_curve : Hazard curve [yearfrac; h]

% separate the coupons
couponDates_1y = couponSchedule_1y(:,1);
couponValues_1y = couponSchedule_1y(:,2);
couponDates_2y = couponSchedule_2y(:,1);
couponValues_2y = couponSchedule_2y(:,2);

% get the zero curve
zeroDates = zeroCurve(:,1);
zeroRates = zeroCurve(:,2);
% interpolate the dates not available in the zero curve
zeroRates_1y = interp1(zeroDates, zeroRates, couponDates_1y);
zeroRates_2y = interp1(zeroDates, zeroRates, couponDates_2y);
% compute discount factors
B_1y = exp(-zeroRates_1y .* couponDates_1y);
B_2y = exp(-zeroRates_2y .* couponDates_2y);

% one year function of the hazard rate
function P_1y = Price_1y(h1)
    survivalProb_1y = [1; exp(-h1 * couponDates_1y)];
    B_bar_1y = B_1y .* exp(-h1 * couponDates_1y);
    P_1y = couponValues_1y' * B_bar_1y + R * 100 * B_1y' * (survivalProb_1y(1:end-1) - survivalProb_1y(2:end)) - dirtyPrice_1y;
end 

P_1y = @(h1) Price_1y(h1);

% compute the hazard rate
h1 = fzero(P_1y, 0.01);

% two year function of the hazard rate
function P_2y = Price_2y(h2, h1)
    % compute the survival probabilities
    survivalProb = zeros(4,1);
    survivalProb(1:2) = exp(-h1 * couponDates_2y(1:2));
    survivalProb(3:4) = exp(-h2 * (couponDates_2y(3:4)-couponDates_2y(2))) * survivalProb(2);
    % compute defaultable bond prices
    B_bar_2y = B_2y .* survivalProb;
    % compute the price
    P_2y = couponValues_2y' * B_bar_2y + R * 100 * B_2y' * ([1; survivalProb(1:end-1)] - survivalProb) - dirtyPrice_2y;
end 

P_2y = @(h2) Price_2y(h2, h1);

h2 = fzero(P_2y, 0.01);

% return the hazard curve in term of basis points
h_curve = [
    couponSchedule_1y(end,1), h1;
    couponSchedule_2y(end,1), h2
];

end