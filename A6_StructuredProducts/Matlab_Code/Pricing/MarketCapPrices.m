function CapPrices = MarketCapPrices(ttms, strikes, Vols, discounts, dates)
% MARKETCAPPRICES computes the cap prices from the market data
%
% INPUTS
%   ttms: time to maturities
%   strikes: strikes
%   Vols: market volatilities
%   discounts: discounts
%   dates: dates of the market data

% save the caps
CapPrices = zeros(size(Vols));

for i = 1:length(ttms)

    % compute the exercise and payment dates
    exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
        calmonths(0:3:12*ttms(i)-3)';
    payment_dates = exercise_dates + calmonths(3);
    % move to business days if needed
    exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
        busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
    payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
        busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
    % convert to datenum
    exercise_dates = datenum(exercise_dates);
    payment_dates = datenum(payment_dates);

    % for each strike
    for j = 1:length(strikes)
        % compute the cap price for given maturity and strike
        CapPrices(i, j) = CapFlat(strikes(j), Vols(i, j), exercise_dates, payment_dates, true, discounts, dates);

    end
    
end

end