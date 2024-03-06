function [ PV ] = PV_risky_bond_h(cf_schedule, h_curve, ZC_curve, R)

l = length(cf_schedule);

if l == 2 % one year case
    % save the coupon dates and values
    couponDates = cf_schedule(:, 1)
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
    PV = couponValues' * B_bar + R * B' * (survivalProb(1:end-1) - survivalProb(2:end));
elseif l == 4 % two year case
    PV = 0;

else
    % throw error
    error('PV_risky_bond_h: invalid cash flow schedule');
end

end