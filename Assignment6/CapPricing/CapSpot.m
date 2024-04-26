function Price = CapSpot(strike, exercise_dates, payment_dates, spot_vols, ttms, strikes, ...
    discounts, dates)
% CapSpot: Compute the price of a cap using the spot volatility assumption
%
% INPUT
%   Strike : Strike price of the cap as a percentage change (eg. 5 for 5%) 
%   exercise_dates : Dates in which the caplets are exercised
%   payment_dates : Dates in which the Euribor is paid
%   spot_vols : Spot volatilities (calibrated to the market caplets)
%   ttms : Time to maturities of the caplets (in years)
%   strikes : Strikes of the caplets
%   discounts : Discount factors curve
%   dates : Dates of the discount factors

% set the date convention
ACT_360 = 2;

% compute the year fractions
delta_ttms = yearfrac(dates(1), exercise_dates, ACT_360);

% initialize the needed spot volatilities
needed_spot_vols = zeros(length(delta_ttms), length(strikes));

% interpolate linearly by columns
for i = 1:length(strikes)
    % interpolate the spot volatilities (extrapolate with flat if needed)
    needed_spot_vols(:,i) = interp1(ttms, spot_vols(:,i), delta_ttms, 'linear', spot_vols(1,i));
end

% initialize the vector of sigmas to use in the Bachelier formula
sigmas = zeros(length(payment_dates), 1);

% interpolate the strike using spline
for i = 1:length(payment_dates)
    sigmas(i) = interp1(strikes, needed_spot_vols(i,:), strike, 'spline');
end 

% compute the discount factors at the payment dates
DF_payment = intExtDF(discounts, dates, payment_dates);

% compute the forward Libors to price the caplets later
DF_exercise = intExtDF(discounts, dates, exercise_dates);
DF_exercise(isnan(DF_exercise)) = 1;
fwd_discounts = DF_payment ./ DF_exercise;
deltas_libor = yearfrac(exercise_dates, payment_dates, ACT_360);
fwd_libor = 1 ./ deltas_libor .* (1 ./ fwd_discounts - 1);

% compute the caplets using the Bachelier formula
Caplets = zeros(length(payment_dates), 1);

for i = 1:length(payment_dates)
    Caplets(i) = CapletBachelier(fwd_libor(i), strike, sigmas(i), payment_dates(i), exercise_dates(i), ...
        dates(1), DF_payment(i));
end

Price = sum(Caplets);

end