function plot_coarse_delta_buckets(buckets, coarse_delta_buckets)
% PLOT_COARSE_DELTA_BUCKETS plots the coarse delta bucket sensitivities
%
% INPUTS
%   buckets: buckets
%   coarse_delta_buckets: coarse delta bucket sensitivities

figure
bar(buckets, coarse_delta_buckets);

% set the labels
xlabel('Buckets');
ylabel('Coarse buckets delta');

% set the title
title('Coarse delta bucket sensitivities');

end