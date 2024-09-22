function FV = FV_risky_bond(IG_cf_schedule_2y, Q, ZC_curve, Recovery) 
% price in one year of the risky bond
%
% INPUTS:
%   IG_cf_schedule_2y: cash flow schedule of the IG bond
%   Q: transition matrix
%   ZC_curve: zero coupon curve
%   Recovery: recovery rate
%
% OUTPUTS:
%   FV: fair value of the bond in one year for IG, HY and Default

FV = zeros(size(Q,1),1);

% interpolate the ZC curve to the cash flow dates
dates = IG_cf_schedule_2y(:,1);
coupons = IG_cf_schedule_2y(:,2);
zeroRates = interp1(ZC_curve(:,1), ZC_curve(:,2), dates, 'linear');
DF = exp(-zeroRates .* dates);

% compute the forward discount factors in one year
future_DF = DF(3:end) / DF(2);

% get the future coupons and dates
future_coupons = coupons(3:end);
future_dates = dates(3:end);

% get the hazard rates from the Q matrix
hazard_rates = -log(1 - Q(:, end));

% IG
Prob_IG = exp( -hazard_rates(1) * (future_dates - dates(2) ) );
B_bar_IG = future_DF .* Prob_IG;
FV(1) = future_coupons' * B_bar_IG + Recovery * 100 * future_DF' * ([1;Prob_IG(1:end-1)] - Prob_IG);

% HY
Prob_HY = exp( -hazard_rates(2) * (future_dates - dates(2) ) );
B_bar_HY = future_DF .* Prob_HY;
FV(2) = future_coupons' * B_bar_HY + Recovery * 100 * future_DF' * ([1;Prob_HY(1:end-1)] - Prob_HY);

% Default
FV(3) =  Recovery * 100;


end