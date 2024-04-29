function sensitivities = vegaBuckets(mkt_vols, ttms, strikes, X_0, spol_A, fixed_rate_B, spol_B, ...
    cap_5y, cap_10y, cap_15y, discounts, dates)
% VEGABUCKETS computes the vega bucket sensitivities for the certificate
%
% INPUTS
%   mkt_vols: market volatilities
%   ttms: time to maturities of the mkt_vols
%   strikes: strikes
%   X_0: upfront payment
%   spol_A: spread over libor for the first quarter
%   fixed_rate_B: fixed rate for the second quarter
%   spol_B: spread over libor for the third quarter
%   cap_5y: strike of the cap from 0 to 5y
%   cap_10y: strike of the cap from 5 to 10y
%   cap_15y: strike of the cap from 10 to 15y
%   discounts: discounts
%   dates: dates of the market data

% initialize the sensitivities
sensitivities = zeros(length(mkt_vols), 1);
% shift is 1 bp
shift = 10^(-4);

% for each maturity, compute the vega bucket sensitivity
for i = 1:length(mkt_vols)
    
    [shift_ttms, shift_vols] = shiftVolsRow(mkt_vols, i, shift, ttms, strikes, discounts, dates);

    % recompute the upfront payment
    X_shift = computeUpfront(shift_vols, shift_ttms, strikes, dates(1), spol_A, fixed_rate_B, spol_B, ...
        cap_5y, cap_10y, cap_15y, discounts, dates);
    % compute the vega bucket sensitivity
    sensitivities(i) = (X_shift - X_0);
end

end