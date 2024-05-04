function spotVols = spotVols(mkt_prices, ttms, strikes, mkt_vols, caplet_ttms, caplet_yf, caplet_DF, Libor)
% spotVols: Compute the spot volatilities of the caps
%
% INPUT
%   mkt_prices : Market prices of the caps
%   ttms : Time to maturities of the caps
%   strikes : Strikes of the caps
%   mkt_vols : Market volatilities of the caps
%   caplet_ttms : Caplet maturities
%   caplet_yf : Caplet year fractions
%   caplet_DF : Caplet discount factors

% initialize the spotVols
spotVols = zeros(length(caplet_ttms), length(strikes));

% first 3 rows is simply the flat vol
spotVols(1:3, :) = repmat(mkt_vols(1, :), 3, 1);

% compute the difference between the cap of following years
Delta_C = mkt_prices(2:end, :) - mkt_prices(1:end-1, :);

for i = 2:length(ttms)

    % find the relevant caplet dates (remember to skip the first)
    relevant_ttms = caplet_ttms(4*ttms(i-1):4*ttms(i)-1);
    relevant_yf = caplet_yf(4*ttms(i-1):4*ttms(i)-1);
    relevant_DF = caplet_DF(4*ttms(i-1):4*ttms(i)-1);
    relevant_Libor = Libor(4*ttms(i-1):4*ttms(i)-1);
    T_alpha = caplet_ttms(4*ttms(i-1)-1);

    for j = 1:length(strikes)

        % find the previous spot vol
        prevVol = spotVols(4*ttms(i-1)-1, j);

        % get the function handle
        fun = @(s) CapSpotBootStrap(strikes(j), prevVol, T_alpha, s, relevant_ttms, relevant_yf, relevant_DF, relevant_Libor) ...
            - Delta_C(i-1, j);

        % compute the spot vol
        sigma_beta = fzero(fun, prevVol);
        % uncomment for faster calibration but less precision
        % sigma_beta = fzero(fun, prevVol, optimset( 'TolX', 1e-6, 'Display', 'off'));

        % compute the spot volatilities
        [~, sigmas] = CapSpotBootStrap(strikes(j), prevVol, T_alpha, sigma_beta, relevant_ttms, relevant_yf, relevant_DF, relevant_Libor);

        % insert into the spotVols
        spotVols(4*ttms(i-1):4*ttms(i)-1, j) = sigmas;
        
    end

end

end