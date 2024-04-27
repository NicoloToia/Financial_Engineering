function plot_hedging_weights(buckets, weights, title_str)
% PLOT_HEDGING_WEIGHTS plots the weights of the hedging
%
% INPUTS
%   buckets: buckets
%   weights: weights of the hedging for each bucket

% create the figure
figure;
% plot the weights vs buckets
bar(buckets, weights);
% set the labels
xlabel('Buckets');
ylabel('Weights');
% set the title
title(title_str);

end