function [alpha, IC] = priceCertificate(S1_0,  sigma_1, d_1, S2_0, sigma_2, d_2, rho, s_spol, P, X, ...
    principal, N, partyA_dates, partyB_dates, dates, discounts, confidence_level)

% priceCertificate: function to compute the participation coefficient of a certificate
%   - party A pays the Euribor 3M + s_spol quarterly with ACT/360 convention and modified following convention
%     and at maturity pays (1-P) * principal
%  - party B pays at start date X * principal and at maturity pays the Coupon. The coupon is
%    computed as the average of the returns of two stocks (ENEL and AXA) each year minus P
%    (if positive) multiplied by alpha.
%    The returns are computed as the ratio of the prices of the stocks.
%
% INPUTS
%   - S1_0: initial price of the first stock
%   - sigma_1: volatility of the first stock
%   - d_1: dividend yield of the first stock
%   - S2_0: initial price of the second stock
%   - sigma_2: volatility of the second stock
%   - d_2: dividend yield of the second stock
%   - rho: correlation between the two stocks
%   - N: number of simulations
%   - fixedLegDates: dates of the fixed leg
%   - floatingLegDates: dates of the floating leg
%   - dates: dates of the discount factors
%   - discounts: discount factors
%
% OUTPUTS
%  - alpha: participation coefficient of the certificate

%% First part: fixed leg (or party B leg)
format long

% compute the year fraction between fixedLegDates
ACT_365 = 3;
deltas_B = yearfrac([dates(1) partyB_dates(1:end-1)], partyB_dates, ACT_365);

% compute the discount factors for the fixed leg
discounts_B = intExtDF(discounts, dates, partyB_dates);
% compute the forward discount factors
forward_discounts_B = discounts_B./ [1 discounts_B(1:end-1)];

% variables to store the simulations
S_1 = zeros(N, length(partyB_dates)+1);
S_2 = zeros(N, length(partyB_dates)+1);
% set initial values
S_1(:,1) = S1_0;
S_2(:,1) = S2_0;

% for each time step
for i = 1:length(partyB_dates)
    % extract the two gaussian random variables
    Z = mvnrnd([0 0], [1 rho; rho 1], N);
    Z_1= Z(:, 1);
    Z_2= Z(:, 2);
    % update the prices
    S_1(:, i+1) = 1/forward_discounts_B(i) * S_1(:, i) .* exp((-d_1- 0.5 * sigma_1^2) * deltas_B(i) + ...
        sigma_1* sqrt(deltas_B(i)) * Z_1);
    S_2(:, i+1) = 1/forward_discounts_B(i) * S_2(:, i) .* exp((-d_2 - 0.5 * sigma_2^2) * deltas_B(i) + ...
        sigma_2 * sqrt(deltas_B(i)) * Z_2);
end 

% compute the returns each year
ratio_1= S_1(:, 2:end) ./ S_1(:, 1:end-1);
ratio_2= S_2(:, 2:end) ./ S_2(:, 1:end-1);

% compute the final coupon
S_t = 1/4 * sum(1/2 * ratio_1 + 1/2 * ratio_2, 2);

% compute the discounted payoffs
disc_payoff = discounts_B(end) * max(S_t - P, 0) * principal;
% compute the mean
maturity_payment_B = mean(disc_payoff);
% initial party B payment
start_payment_B = X * principal;

%% Party A leg (or floating leg)

% compute discount factors for the floating leg
discounts_A = intExtDF(discounts, dates, partyA_dates);

% compute the deltas
ACT_360 = 2;
deltas_A = yearfrac([dates(1) partyA_dates(1:end-1)], partyA_dates, ACT_360);

% compute the s_spol payments
NPV_s_spol = s_spol * sum(deltas_A.* discounts_A) * principal;

% compute the Libor payments
NPV_libor = (1 - discounts_A(end)) * principal;

% payment at maturity
maturity_payment_A = (1-P) * discounts_A(end) * principal;

party_A_leg = NPV_libor + NPV_s_spol + maturity_payment_A;

%% Compute the partecipation coefficient

alpha = (party_A_leg - start_payment_B) / maturity_payment_B;

%% Compute the confidence interval

% compute the standard deviation
std_alpha = std((NPV_libor + NPV_s_spol + disc_payoff - start_payment_B) / maturity_payment_B);
% compute the confidence interval
z = norminv(1 - (1 - confidence_level) / 2);
IC = [alpha - z * std_alpha / sqrt(N), alpha + z * std_alpha / sqrt(N)];

end