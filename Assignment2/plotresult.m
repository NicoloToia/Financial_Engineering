function plotresult(dates, discounts, zeroRates)
% plotresult: Plots the discount factors and zero rates
%
% INPUTS:
% dates: dates of the discount factors and zero rates
% discounts: discount factors
% zeroRates: zero rates
%

% create a new figure
figure 

% left hand side
yyaxis left
% plot discount factors as green filled triangles
plot(dates(1:end), discounts(1:end), 'g-^', 'MarkerFaceColor', 'g')
ylabel('Discount Factors')
ylim([0 1]) % make the y-axis go from 0 to 1
yticks(0:0.2:1) % ticks every 0.2

%zero-rates
yyaxis right
% plot zero rates as blue filled diamonds (in percentage)
plot(dates(1:end), zeroRates(1:end), 'b-d', 'MarkerFaceColor', 'b')
ylim([2.5 5.0]) % make the y-axis go from 2.5 to 5.0
yticks(2.5:0.5:5.0) % ticks every 0.5
ylabel('Zero Rates')

% set x-axis limits
xlabel('Date')
title('Zero Rates and Discount Factors')
legend('Discount Factors', 'Zero Rates')
grid on
% x-axis in date format Month name and last two digits of year
datetick('x', 'mmm-yy')
% set the x-axis limits (one year before and after the dates)
xlim([dates(1)-365 dates(end)+365])
% make axes black
ax = gca;
ax.XColor = 'k';
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

end