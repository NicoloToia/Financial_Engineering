function PlotIntensities(datesCDS, int_Approx, int_Exact, int_JT)
% Plot of the Bootstrap results considering all the different flag case:
% 1 approx; 2 exact; 3 Jarrow-Turnbull
%
% INPUT
%   datesCDS    : dates of the CDS
%   int_Approx  : intesities neglecting the accrual
%   int_Exact   : intensities considering the accrual
%   int_JT      : flat intensities Jarrow-Turmbull
%
% OUTPUT
%   Plots of int_Approx vs int_Exact & qualitative plot int_JT

% Plot the intensities with and without accrual
figure
stairs(datesCDS, int_Approx, '.-')
hold on
stairs(datesCDS, int_Exact, '.-')
legend('Approx', 'Exact')
datetick('x', 'mm/dd/yyyy', 'keepticks')

% Compute the cumulative mean of the intensities
cumIntensities = cumsum(int_Approx) ./ (1:length(int_Approx))';

% Plot of cumulative mean and JT intensities (qualitative), int_JT are flat
figure
plot(datesCDS, cumIntensities, 'o-')
hold on
stairs(datesCDS, int_Approx, '.-')
stairs(datesCDS, int_JT, '-')
legend('Mean', 'Approx', 'JT')
grid on
xlabel('Date')
ylabel('Intensity')
datetick('x', 'mm/dd/yyyy', 'keepticks')

end