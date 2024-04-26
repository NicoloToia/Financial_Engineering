function plot_vega_buckets(vega_buckets, ttms, strikes)

% create the figure
figure;
% plot the vega_bucket sensitivities as a surface
surf(strikes, ttms, vega_buckets);
% set the labels
xlabel('Time to maturity (years)');
ylabel('Strikes');
zlabel('Vega bucket sensitivity');
% set the title
title('Vega bucket sensitivities');
% set the grid
grid on;

% set the colorbar
colorbar;

end
