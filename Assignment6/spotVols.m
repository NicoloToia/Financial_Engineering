function spotVols = spotVols(mkt_prices, ttms, strikes, mkt_vols, discounts, dates)
% spotVols: Compute the spot volatilities of the caps
%
% INPUT
%   mkt_prices : Market prices of the caps
%   ttms : Time to maturities of the caps
%   strikes : Strikes of the caps
%   mkt_vols : Market volatilities of the caps
%   discounts : Discount factors
%   dates : Dates of the discount factors

% initialize the spotVols
spotVols = zeros(length(ttms), length(strikes));

% first row is simple the flat vol
spotVols(1, :) = mkt_vols(1, :);

% compute the difference between the cap of following years
diffCap = mkt_prices(2:end, :) - mkt_prices(1:end-1, :);

for i = 2:length(ttms)

    % compute the start date
    start_date = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
        calyears(1*ttms(i-1))';
    % compute the exercise dates
    exercise_dates = start_date + calmonths(0:3:12*(ttms(i)-ttms(i-1))-3)';
    exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
        busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
    % compute the payment dates
    payment_dates = start_date + calmonths(3:3:12*(ttms(i)-ttms(i-1)))';
    payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
        busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());

    for j = 1:length(strikes)

        % get the function handle
        fun = @(s) CapSpotBootStrap(strikes(j), spotVols(i-1, j), s, startDates(i-1), paymentDates, discounts, dates) - ...
            diffCap(i-1, j);

        % compute the spot vol with fmincon (with positive constraint)
        spotVols(i, j) = fzero(fun, spotVols(i-1, j));
        
    end
end

end