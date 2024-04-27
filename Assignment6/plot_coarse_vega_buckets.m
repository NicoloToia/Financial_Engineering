function plot_coarse_vega_buckets(buckets, coarse_vega_buckets)
% PLOT_COARSE_VEGA_BUCKETS plots the coarse vega bucket sensitivities
%
% INPUTS
%   buckets: buckets
%   coarse_vega_buckets: coarse vega bucket sensitivities

figure
bar(buckets, coarse_vega_buckets);

% set the labels
xlabel('Buckets');
ylabel('Coarse buckets vega');

% set the title
title('Coarse vega bucket sensitivities');

end