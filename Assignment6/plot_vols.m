function plotVols(spot_vols, spot_ttms, mkt_vols, ttms, strikes)
% plot the spot volatilities surface against the flat volatilities
%
% INPUTS
%   spot_vols: spot volatilities
%   mkt_vols: market volatilities
%   ttms: time to maturities
%   strikes: strikes

figure;

% plot the spot volatilities surface against the flat volatilities
% spot vols in blue, market vols in red
surf(strikes, spot_ttms, spot_vols, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'blue');
hold on
surf(strikes, ttms, mkt_vols, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'FaceColor', 'red');
xlabel('Strikes');
ylabel('TTMs');
zlabel('Volatilities');
title('Spot Volatilities vs Market Volatilities');
legend('Spot Volatilities', 'Market Volatilities');

end