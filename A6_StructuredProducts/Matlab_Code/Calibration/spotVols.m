function [spot_ttms, spotVols] = spotVols(mkt_prices, ttms, strikes, mkt_vols, discounts, dates)
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
spotVols = zeros(4*ttms(end)-1, length(strikes));

% first 3 rows is simply the flat vol
spotVols(1:3, :) = repmat(mkt_vols(1, :), 3, 1);

% compute the difference between the cap of following years
diffCap = mkt_prices(2:end, :) - mkt_prices(1:end-1, :);

% start the waitbar
wb = waitbar(0, 'Computing the spot volatilities...');
total = length(ttms);

for i = 2:length(ttms)

    % compute the start date
    start_date = datetime(dates(1), 'ConvertFrom', 'datenum') + ...
        calyears(ttms(i-1))';
    % compute the T alpha (the last exercise date of the previous cap)
    T_alpha = start_date - calmonths(3)';
    % compute the exercise dates and payment dates
    exercise_dates = start_date + calmonths(0:3:12*(ttms(i)-ttms(i-1))-3)';
    payment_dates = exercise_dates + calmonths(3)';
    % move to business days
    if ~isbusday(T_alpha, eurCalendar())
        T_alpha = busdate(T_alpha, 'modifiedfollow', eurCalendar());
    end
    exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
        busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
    payment_dates(~isbusday(payment_dates, eurCalendar())) = ...
        busdate(payment_dates(~isbusday(payment_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
    % convert to datenum
    T_alpha = datenum(T_alpha);
    exercise_dates = datenum(exercise_dates);
    payment_dates = datenum(payment_dates);

    for j = 1:length(strikes)

        % find the previous spot vol
        prevVol = spotVols(4*ttms(i-1)-1, j);

        % get the function handle
        fun = @(s) CapSpotBootStrap(strikes(j), prevVol, T_alpha, s, exercise_dates, payment_dates, discounts, dates) - ...
            diffCap(i-1, j);

        % compute the spot vol
        sigma_beta = fzero(fun, prevVol);
        % uncomment for faster calibration but less precision
        % sigma_beta = fzero(fun, prevVol, optimset( 'TolX', 1e-6, 'Display', 'off'));

        % compute the spot volatilities
        [~, sigmas] = CapSpotBootStrap(strikes(j), prevVol, T_alpha, sigma_beta, exercise_dates, payment_dates, discounts, dates);

        % insert into the spotVols
        spotVols(4*ttms(i-1):4*ttms(i)-1, j) = sigmas;
        
    end

    % update the waitbar with the percentage
    waitbar(i/total, wb, sprintf('Computing the spot volatilities... %.2f%%', i/total*100));
end

% close the waitbar
close(wb);

% compute the new time to maturities
exercise_dates = datetime(dates(1), 'ConvertFrom', 'datenum') + calmonths(3:3:12*ttms(end)-3)';
% move to business days
exercise_dates(~isbusday(exercise_dates, eurCalendar())) = ...
    busdate(exercise_dates(~isbusday(exercise_dates, eurCalendar())), 'modifiedfollow', eurCalendar());
exercise_dates = datenum(exercise_dates);
ACT_365 = 3;
spot_ttms = yearfrac(dates(1), exercise_dates, ACT_365);

end