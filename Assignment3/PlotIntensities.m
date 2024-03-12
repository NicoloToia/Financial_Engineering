function PlotIntensities(datesCDS, int_Approx, int_Exact, int_JT)

    figure
    stairs(datesCDS, int_Approx, '-')
    hold on
    stairs(datesCDS, int_Exact, '-')
    legend('Approx', 'Exact')
    
    
    % compute the cumulative mean of the intensities
    cumIntensities = cumsum(int_Approx) ./ (1:length(int_Approx))';
    % plot the intensities as a step function
    figure
    plot(datesCDS, cumIntensities, '-')
    hold on
    stairs(datesCDS, int_Approx, '-')
    stairs(datesCDS, int_JT, '-')
    legend('Mean', 'Approx', 'JT')
    grid on
    xlabel('Date')
    ylabel('Intensity')
    datetick('x', 'mm/dd/yyyy', 'keepticks')
end