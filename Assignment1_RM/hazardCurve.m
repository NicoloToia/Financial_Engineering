function [h_curve] = hazardCurve(zeroCurve, R, dirtyPrice_1y, dirtyPrice_2y, ...
    couponSchedule_1y, couponSchedule_2y)

% get the zero curve
zeroDates = zeroCurve(:,1);
zeroRates = zeroCurve(:,2);
% interpolate the dates not available in the zero curve
zeroRates = interp1(zeroDates, zeroRates, couponSchedule_2y(:,1));
discountFactors = exp(-zeroRates .* couponSchedule_2y(:,1));

% one year function of the hazard rate
coupons_1y = couponSchedule_1y(:,2);
yearfrac_1y = [0; couponSchedule_1y(:,1)];
discountFactors_1y = discountFactors(1:2);
discountCoupon_1y = discountFactors_1y .* coupons_1y;
% price as a function of the hazard rate
P_1y = @(h1) exp(-h1 * yearfrac_1y(2:3))' * discountCoupon_1y + ...
    R * ( exp(-h1 * yearfrac_1y(1:end-1)) - exp(-h1 * yearfrac_1y(2:end)) )' * discountFactors_1y - ...
    dirtyPrice_1y;
% compute the hazard rate
h1 = fzero(P_1y, 0.01);

% two year function of the hazard rate
coupons_2y = couponSchedule_2y(:,2);
yearfrac_2y = [0; couponSchedule_2y(:,1)];
discountFactors_2y = discountFactors(1:4);
discountCoupon_2y = discountFactors_2y .* coupons_2y;
prevProb = exp(-h1 * yearfrac_2y(3));
% price as a function of the hazard rate
P_2y = @(h2) exp(-h1 * yearfrac_2y(2:3))' * discountCoupon_2y(1:2) + ...
    prevProb * exp(-h2 * (yearfrac_2y(4:5) - yearfrac_2y(3)))' *  discountCoupon_2y(3:4) + ...
    R * ( ...
        ( exp(-h1 * yearfrac_2y(1:2)) - exp(-h1 * yearfrac_2y(2:3)) )' * discountFactors_2y(1:2) + ...
        prevProb * (exp(-h2 * yearfrac_2y(4:5)) - exp(-h2 * yearfrac_2y(3)))' ...
        * discountFactors_2y(3:4) ...
    ) - dirtyPrice_2y;

h2 = fzero(P_2y, 0.01);

% return the hazard curve in term of basis points
h_curve = [
    yearfrac_1y(3), h1;
    yearfrac_2y(5), h2;
];

end