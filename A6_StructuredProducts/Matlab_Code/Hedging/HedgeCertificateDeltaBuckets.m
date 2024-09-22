function weights = HedgeCertificateDeltaBuckets(buckets, delta_buckets, delta_buckets_swaps, doplot)
% HEDGECERTIFICATEDELTABUCKETS hedging of the delta of the certificate
%
% INPUTS
%   buckets: buckets
%   delta_buckets: delta bucket sensitivities
%   delta_buckets_swaps: delta bucket sensitivities for the swaps
%   doplot: plot the hedging steps

% initialize the weights and the portfolio delta
portfolio_delta = delta_buckets;
weights = zeros(size(buckets));

% open the figure if needed
if doplot
    figure
end

% for each bucket
for i = length(buckets):-1:1

    % find the weight to perfectly hedge the bucket using the corresponding swap
    weights(i) = - portfolio_delta(i) / delta_buckets_swaps(i,i);

    % update the portfolio delta
    portfolio_delta = portfolio_delta + weights(i) * delta_buckets_swaps(:,i);

    % subplot the steps of the portfolio delta
    if doplot
        subplot(2,2, length(buckets)-i+1);
        plot(buckets, portfolio_delta*50*10^6, 'o-', 'LineWidth', 2);
        hold on
        % zero line
        plot([0, 15], [0, 0], 'k--', 'LineWidth', 1);
        title(['Portfolio Delta hedged with ', num2str(buckets(i)), ' years swap']);
        xlabel('Years');
        ylabel('Delta');
        % fix the y-axis scale
        ylim([-50000, 15000]);
        grid on;
    end

end

end