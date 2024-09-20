function plot_Vega_matrix(vega_matrix, ttms, strikes)
% PLOT_VEGA_MATRIX plots the vega matrix
%
% INPUTS
%   vega_matrix: vega matrix
%   ttms: time to maturities
%   strikes: strikes


figure;
surf(strikes, ttms, vega_matrix)
xlabel('Strikes')
ylabel('ttms')
zlabel('Vega matrix')
title('Vega matrix')
colormap('parula') % Change the colormap to 'parula' for better visibility



end