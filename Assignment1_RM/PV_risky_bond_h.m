function [ PV ] = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R)
% PV_risky_bond_h: computes the present value of a risky bond using the
% survival probabilities (derived from the hazard rates) and the zero
% coupon curve
%
% INPUT
%   cf_schedule: a 2xN matrix with the cash flow schedule of the bond
%   h_curve: a 2xN matrix with the hazard rates
%   ZC_curve: a 2xN matrix with the zero coupon curve
%   R: recovery rate
%
% OUTPUT
%   PV : Dirty price for a given risky bond (from hazard rate curve - fixed recoovery rate R)

% check the length of the cash flow schedule
l = length(cf_schedule);

if l == 2 % one year case
    % save the coupon dates and values
    couponDates = cf_schedule(:, 1);
    couponValues = cf_schedule(:, 2);
    % compute the discount factors
    ZC_dates = ZC_curve(:, 1);
    ZC_rates = ZC_curve(:, 2);
    ZC_curve = interp1(ZC_dates, ZC_rates, couponDates);
    B = exp(-ZC_curve .* couponDates);
    % compute the survival probabilities
    h_curve = h_curve(1, 2);
    survivalProb = exp(-h_curve * couponDates);
    B_bar = B .* survivalProb;
    % add the survival probability of 1 at time 0
    survivalProb = [1; survivalProb];
    PV = couponValues' * B_bar + R * 100 * B' * (survivalProb(1:end-1) - survivalProb(2:end));
elseif l == 4 % two year case
    % save the coupon dates and values
    couponDates = cf_schedule(:, 1);
    couponValues = cf_schedule(:, 2);
    % compute the discount factors
    ZC_dates = ZC_curve(:, 1);
    ZC_rates = ZC_curve(:, 2);
    ZC_curve = interp1(ZC_dates, ZC_rates, couponDates);
    % compute the survival probabilities
    h_curve = h_curve(:, 2);
    survivalProb = zeros(4, 1);
    survivalProb(1:2) = exp(-h_curve(1) * couponDates(1:2));
    survivalProb(3:4) = exp(-h_curve(1) * couponDates(2)) * ...
        exp(-h_curve(2) * (couponDates(3:4) - couponDates(2)));
    B = exp(-ZC_curve .* couponDates);
    B_bar = B .* survivalProb;
    % add the survival probability of 1 at time 0
    survivalProb = [1; survivalProb];
    PV = couponValues' * B_bar + R * 100 * B' * (survivalProb(1:end-1) - survivalProb(2:end));
else
    % throw error
    error('PV_risky_bond_h: invalid cash flow schedule');
end

end