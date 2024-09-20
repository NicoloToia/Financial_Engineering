function [ttms, strikes, Vols] = readVolData(filename)
% Reads data from excel
%  It reads flat volatilities for different maturities and strikes
%  All volatilities are in bp units
%
% INPUTS:
%  filename: excel file name where data are stored
%  formatData: data format in Excel
% 
% OUTPUTS:
%  dates: maturities
%  strikes: strikes
%  Vols: volatilities

%% Dates from Excel

% read the needed cells
ttms = readtable(filename, 'Sheet', 1, 'Range', 'B2:B17', 'ReadVariableNames', false);
% convert to array
ttms = table2array(ttms);
% remove the second element
ttms = [ttms(1); ttms(3:end)];
% convert from string to number (delete last character (y))
ttms = str2double(cellfun(@(x) x(1:end-1), ttms, 'UniformOutput', false));

%% Strikes from Excel

strikes = readtable(filename, 'Sheet', 1, 'Range', 'F1:R1', 'ReadVariableNames', false);
strikes = table2array(strikes);

%% Vols from Excel

Vols = readtable(filename, 'Sheet', 1, 'Range', 'F2:R17', 'ReadVariableNames', false);
Vols = table2array(Vols);
Vols = Vols / 10000;
% delete the second row
Vols = [Vols(1,:); Vols(3:end,:)];

end % readVolData