function CapPrices = MarketCapPrices(ttms, strikes, Vols, discounts, dates)
%

% save the caps
CapPrices = zeros(length(ttms), length(strikes));

for i = 1:length(ttms)
    % compute the exercise dates
    exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
        calmonths(0:3:12*ttms(i)-3)';
    % move to business days if needed
    exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
        busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar())
    exercise_dates = datenum(exercise_dates);

    % compute the payment dates
    payment_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
        calmonths(3:3:12*ttms(i))';
    payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
        busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar())
    payment_dates = datenum(payment_dates);

    for j = 1:length(strikes)

        CapPrices(i, j) = CapFlat(strikes(j), Vols(i, j), exercise_dates, payment_dates, true, discounts, dates);

    end
    
end

end