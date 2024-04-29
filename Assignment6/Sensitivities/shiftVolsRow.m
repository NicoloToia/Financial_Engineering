function [shift_ttms, shift_vols] = shiftVolsRow(mkt_vols, target_row, shift, ttms, strikes, discounts, dates)
% SHIFTVOLSROW shifts the row of the market volatilities by a given shift
%
% INPUTS
%   mkt_vols: market volatilities
%   target_row: row to shift
%   shift: shift
% OUTPUTS
%   shift_ttms: time to maturities of the shifted volatilities
%   shift_vols: shifted spot volatilities
%

% check if the shift has already been calibrated for in the file
if isfile('Data/calibrated_vols.mat')
    load('Data/calibrated_vols.mat', 'calibrated_vols')
else
    % create the shifted vols cell array
    calibrated_vols = cell(length(mkt_vols), 1)
    save('Data/calibrated_vols.mat', 'calibrated_vols')
end

% check that the row has not been calibrated for yet and fill the cell array
if isempty(calibrated_vols{target_row})
    % shift the row by 1 bp
    shift_mkt_vols = mkt_vols;
    shift_mkt_vols(target_row, :) = shift_mkt_vols(target_row, :) + shift;
    % recompute the cap prices
    mkt_prices = MarketCapPrices(ttms, strikes, shift_mkt_vols, discounts, dates);
    % recalibrate the spot vols
    [shift_ttms, shift_vols] = spotVols(mkt_prices, ttms, strikes, shift_mkt_vols, discounts, dates);

    % save the shifted vols
    calibrated_vols{target_row} = {shift_ttms, shift_vols};

    % save the calibrated vols
    save('Data/calibrated_vols.mat', 'calibrated_vols');
end

% retrieve the shifted vols
shift_ttms = calibrated_vols{target_row}{1};
shift_vols = calibrated_vols{target_row}{2};

end