function plot_coarse_delta_buckets_swaps(buckets, coarse_delta_buckets_swaps)
% PLOT_COARSE_DELTA_BUCKETS_SWAPS plots the coarse delta bucket sensitivities for the swaps
%
% INPUTS
%   buckets: buckets
%   coarse_delta_buckets_swaps: coarse delta bucket sensitivities for the swaps

figure
subplot(2,2,1)
bar(buckets, coarse_delta_buckets_swaps(:,1));
% set the labels
xlabel('Buckets');
ylabel('Coarse buckets delta');
% set the title
title('Coarse delta bucket sensitivities for 2 years swap');

subplot(2,2,2)
bar(buckets, coarse_delta_buckets_swaps(:,2));
% set the labels
xlabel('Buckets');
ylabel('Coarse buckets delta');
% set the title
title('Coarse delta bucket sensitivities for 5 years swap');

subplot(2,2,3)
bar(buckets, coarse_delta_buckets_swaps(:,3));
% set the labels
xlabel('Buckets');
ylabel('Coarse buckets delta');
% set the title
title('Coarse delta bucket sensitivities for 10 years swap');

subplot(2,2,4)
bar(buckets, coarse_delta_buckets_swaps(:,4));
% set the labels
xlabel('Buckets');
ylabel('Coarse buckets delta');
% set the title
title('Coarse delta bucket sensitivities for 15 years swap');

end