function plot_vega_buckets(vega_buckets, ttms)
% PLOT_VEGA_BUCKETS plots the vega bucket sensitivities
%
% INPUTS
%   vega_buckets: vega bucket sensitivities
%   ttms: time to maturities

% create the figure
figure;
% plot the vega_bucket sensitivities vs time to maturity
plot(ttms, vega_buckets, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5);
% set the labels
xlabel('Time to maturity (years)');
ylabel('Vega bucket sensitivity');
% set the title
title('Vega bucket sensitivities');
% set the grid
grid on;

end