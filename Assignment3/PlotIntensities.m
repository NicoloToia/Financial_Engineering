function PlotIntensities(datesCDS, intensities, t0)

figure
dates=[t0;datesCDS];
for i=1:length(intensities)
    xi=linspace(dates(i),dates(i+1));
    yi=ones(length(xi),1)*intensities(i);
    hold on
    plot(xi,yi,"LineWidth",3,"Color","b")
end
hold off