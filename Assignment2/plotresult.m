function plotresult(dates, discounts, zeroRates)

%discount curve
figure 
yyaxis left
% plot discount factors as green filled triangles
plot(dates(2:end), discounts(2:end), 'g-^', 'MarkerFaceColor', 'g')
ylabel('Discount Factors')
% make the y-axis go from 0 to 1
ylim([0 1])
% ticks every 0.2
yticks(0:0.2:1)

%zero-rates
yyaxis right
% plot zero rates as blue filled diamonds (in percent)
plot(dates(2:end), zeroRates(2:end), 'b-d', 'MarkerFaceColor', 'b')
% make the y-axis go from 2.5 to 5.0
ylim([2.5 5.0])
yticks(2.5:0.5:5.0)
ylabel('Zero Rates')

% set x-axis label and title
xlabel('Date')
title('Zero Rates and Discount Factors')
legend('Discount Factors', 'Zero Rates')
grid on
% x-axis in date format Month name and last two digits of year
datetick('x', 'mmm-yy')
% cut off dates and rates after the Jan-2042
% make axes black
ax = gca;
ax.XColor = 'k';
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';