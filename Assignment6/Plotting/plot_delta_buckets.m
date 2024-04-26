function plot_delta_buckets(dates, buckets)
% PLOT_DELTA_BUCKETS plots the delta-bucket sensitivities
%
% INPUTS
%   dates: dates of the market data
%   buckets: delta-bucket sensitivities

figure;
plot(dates, buckets, 'LineWidth', 2);
title('Delta-Bucket Sensitivities');
xlabel('Dates');
ylabel('Delta-Bucket Sensitivities');
datetick('x', 'mm/dd/yyyy', 'keeplimits', 'keepticks');
grid on;

end