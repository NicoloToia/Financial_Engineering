function [dates, rates] = readRatesData( filename, formatData)
% Reads data from excel
%  It reads bid/ask prices and relevant dates
%  All input rates are in % units
%
% INPUTS:
%  filename: excel file name where data are stored
%  formatData: data format in Excel
% 
% OUTPUTS:
%  dates: struct with settlementDate, deposDates, futuresDates, swapDates
%  rates: struct with deposRates, futuresRates, swapRates

%% Dates from Excel

%Settlement date
[~, settlement] = xlsread(filename, 1, 'E8');
%Date conversion
dates.settlement = datenum(settlement, formatData);

%Dates relative to depos
[~, date_depositi] = xlsread(filename, 1, 'D11:D21');
dates.depos = datenum(date_depositi, formatData);

%Dates relative to futures: calc start & end
[~, date_futures_read] = xlsread(filename, 1, 'Q12:R19');
numberFutures = size(date_futures_read,1);

dates.futures=ones(numberFutures,2);
dates.futures(:,1) = datenum(date_futures_read(:,1), formatData);
dates.futures(:,2) = datenum(date_futures_read(:,2), formatData);

% initialize swaps as empty
dates.swaps = [];

%% Rates from Excel (Bids & Asks)

%Depos
tassi_depositi = xlsread(filename, 1, 'E11:F21');
rates.depos = tassi_depositi / 100;
% get the mid market rate
rates.depos = mean(rates.depos, 2);

%Futures
tassi_futures = xlsread(filename, 1, 'E28:F35');
%Rates from futures
tassi_futures = 100 - tassi_futures;
rates.futures = tassi_futures / 100;
% get the mid market rate
rates.futures = mean(rates.futures, 2);

%Swaps
tassi_swaps = xlsread(filename, 1, 'E39:F56');
rates.swaps = tassi_swaps / 100;
% get the mid market rate
rates.swaps = mean(rates.swaps, 2);

end % readExcelData