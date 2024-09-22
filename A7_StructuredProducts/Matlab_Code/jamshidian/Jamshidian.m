function price = Jamshidian(alpha, ttm, strike_swaption, a, sigma, discounts, dates)
%   Jamshidian: Computes the price of a swaption using the Jamshidian formula
%
% INPUTS
%   alfa: year fraction of the first exercise date
%   omega: year fraction of the expiry date
%   strike_swaption: strike of the swaption
%   a: mean reversion of the Hull-White model
%   sigma: volatility of the Hull-White model
%   discounts: discount factors
%   dates: dates for the discount factors

% define the year fraction conventions
ACT_360 = 2;
ACT_365 = 3;

% compute the start date and the dates in which the bonds are paid
T_alpha = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(alpha)';
if ~isbusday(T_alpha, eurCalendar())
    T_alpha = busdate(T_alpha, 'modifiedfollow', eurCalendar());
end
T_alpha = datenum(T_alpha);

bond_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(alpha+1:ttm)';
bond_dates(~isbusday(bond_dates, eurCalendar())) = ...
    busdate(bond_dates(~isbusday(bond_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
bond_dates = datenum(bond_dates);

% compute the forward discount factors from t_alpha to each payment date (evaluated at t0)
B_alpha = intExtDF(discounts, dates, T_alpha);
fwd_ZCB_t0 = intExtDF(discounts, dates, bond_dates) / B_alpha;

% find the x that satisfies the coupon bond equation
x_star = fzero(@(x) couponBond(x, fwd_ZCB_t0, a, sigma, dates(1), T_alpha, bond_dates, strike_swaption) - 1, 0);

% compute the strikes for each zcb call using the x_star
K = zeros(length(bond_dates), 1);
ti = yearfrac(dates(1), T_alpha, ACT_365);
for i = 1:length(bond_dates)
    % compute the year fraction
    tau = yearfrac(T_alpha, bond_dates(i), ACT_365);
    K(i) = fwd_ZCB_t0(i) * exp(-x_star * sigmaHJM(a, sigma, tau, 0)/sigma - 0.5 * IntHJM(a, sigma, ti, 0, tau));
end


% compute the coupons
coupon_deltas = yearfrac([T_alpha; bond_dates(1:end-1)], bond_dates, ACT_360);
coupons = strike_swaption * coupon_deltas;
coupons(end) = coupons(end) + 1;

% for each K_i compute the ZCB call option
ZCB_call = zeros(length(bond_dates), 1);
for i = 1:length(bond_dates)
    ZCB_call = ZC_call(K(i), fwd_ZCB_t0(i), B_alpha, T_alpha, dates(1), sigma, a, bond_dates(i));
end

% coupon call
coupon_call = sum(coupons .* ZCB_call);

% coupon bond 
coupon_bond = sum(coupons .* fwd_ZCB_t0) * B_alpha;

% put-call parity
price = coupon_call + B_alpha * 1 - coupon_bond;

end