function plot_full_vega_buckets(vega_buckets, ttms, strikes)
% PLOT_VEGA_BUCKETS plots the vega bucket sensitivities surface
%
% INPUTS
%   vega_buckets: vega bucket sensitivities
%   ttms: time to maturities
%   strikes: strikes

% print the vega buckets sensitivities as a table

% create the row names from the ttms
row_names = cell(1, length(ttms));
for i = 1:length(ttms)
    row_names{i} = ['T = ', num2str(ttms(i)), 'y'];
end

% create the column names from the strikes
column_names = cell(1, length(strikes));
for i = 1:length(strikes)
    column_names{i} = ['K = ', num2str(strikes(i))];
end

% create the table
vega_buckets_table = array2table(vega_buckets, 'RowNames', row_names, 'VariableNames', column_names);

% print the table
disp(vega_buckets_table);

% create the figure
figure;
% plot the vega_bucket sensitivities as a surface
% TODO: invert the names of y and x axis
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
