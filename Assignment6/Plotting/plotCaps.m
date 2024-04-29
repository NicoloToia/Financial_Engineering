function plotCaps(Caps, ttms, strikes)
    % PLOTCAPS plots the cap surface
    %
    % INPUTS
    %   Caps: cap prices
    %   ttms: time to maturities
    %   strikes: strikes
    
    figure;

    % plot the spot volatilities surface against the flat volatilities
    % spot vols in blue, market vols in red
    surf(strikes, ttms, Caps, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'blue');
    xlabel('Strikes');
    ylabel('TTMs');
    zlabel('Cap Price');
    title('Cap Surface');
    
    end