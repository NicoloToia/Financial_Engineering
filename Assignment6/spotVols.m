function spotVols = spotVols(ttms, strikes, capPrices, Vols, discounts, dates)

% initialize the spotVols
spotVols = zeros(length(ttms), length(strikes));

% first row is simple the flat vol
spotVols(1, :) = Vols(1, :);

% compute the difference between the cap of following years
diffCap = capPrices(2:end, :) - capPrices(1:end-1, :);

% compute the starting date for each cap
startDates = datetime(dates(1), 'ConvertFrom', 'datenum') + calyears(ttms)';
% move to business days if needed
startDates(~isbusday(startDates, eurCalendar())) = ...
    busdate(startDates(~isbusday(startDates, eurCalendar())), 'modifiedfollow', eurCalendar())
startDates = datenum(startDates);

for i = 2:length(ttms)

    % compute the paymentDates
    paymentDates = datetime(startDates(i-1), 'ConvertFrom', 'datenum') + ...
        calmonths(3:3:12*(ttms(i) - ttms(i-1)))';
    paymentDates(~isbusday(paymentDates, eurCalendar())) = ...
        busdate(paymentDates(~isbusday(paymentDates, eurCalendar())), 'modifiedfollow', eurCalendar())
    paymentDates = datenum(paymentDates);

    for j = 1:length(strikes)

        % get the function handle
        fun = @(s) CapSpot(strikes(j), spotVols(i-1, j), s, startDates(i-1), paymentDates, false, discounts, dates) - diffCap(i-1, j);

        % compute the spot vol with fzero
        spotVols(i, j) = fzero(fun, Vols(i, j));
        
    end
end


end