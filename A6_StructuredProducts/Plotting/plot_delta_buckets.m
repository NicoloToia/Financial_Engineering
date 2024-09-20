function plot_delta_buckets(dates, buckets)
% PLOT_DELTA_BUCKETS plots the delta-bucket sensitivities
%
% INPUTS
%   dates: dates of the market data
%   buckets: delta-bucket sensitivities

figure;
% plot the points (filled circles) and connect them with a line
plot(dates, buckets, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5);
title('Delta-Bucket Sensitivities');
xlabel('Dates');
ylabel('Delta-Bucket Sensitivities');
% show all dates rotated by 45 degrees
datetick('x', 'mm/dd/yyyy', 'keeplimits', 'keepticks');
% date_ticks = datestr(dates, 'dd-mmm-yyyy');
% set(gca, 'XTick', dates, 'XTickLabel', date_ticks, 'XTickLabelRotation', 45);

grid on;

end