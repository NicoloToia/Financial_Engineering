function price = couponBond(x, fwd_DF_t0, a, sigma, t0, t_alpha, bond_dates, strike_swaption)
%  couponBond: Computes the price of a coupon bond using the Hull-White model and Jamshidian formula
%
% INPUTS
%   x: the stochastic interest rate
%   fwd_DF_t0: forward discount factor from the start date to the exercise date evaluated at t0
%   a: mean reversion of the Hull-White model
%   sigma: volatility of the Hull-White model
%   t0: settlement date of the swaption
%   t_alpha: start date of the bond
%   bond_dates: dates in which the coupons are paid
%   strike_swaption: strike of the swaption

% define the year fraction conventions
ACT_360 = 2;

% compute the coupons and dates
coupon_deltas = yearfrac([t_alpha; bond_dates(1:end-1)], bond_dates, ACT_360);
coupons = strike_swaption * coupon_deltas;
coupons(end) = coupons(end) + 1;

% compute the forward discounts based on the interest rate x
fwd_DF = zeros(length(bond_dates), 1);
% extreme of integration for the sigma
ti = yearfrac(t0, t_alpha, ACT_360);

% compute the fwd discount factors
for i = 1:length(fwd_DF)
    % compute the relevant year fractions
    tau = yearfrac(t_alpha, bond_dates(i), ACT_360);
    % compute the forward discount factor using lemma 1
    fwd_DF(i) = fwd_DF_t0(i) * exp(-x * sigmaHJM(a, sigma, tau, 0)/sigma - 0.5 * IntHJM(a, sigma, ti, 0, tau));
end

% compute the price of the bond
price = sum(coupons .* fwd_DF);

end