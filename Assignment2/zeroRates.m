function [zRates] = zeroRates(dates, discounts)
% zeroRates: Compute zero rates from discount factors computed in function bootstrap
%
% INPUT
% dates    : Dates of the discount factors
% discounts: Discount factors computed in the bootstrap
%
% OUTPUT:
% zRates   : Zero rates

% Compute zero rates 
zRates = -log(discounts) ./ yearfrac(dates(1), dates);
% Set the first zero rate equal to the first depo (O/N)
zRates(1) = zRates(2);
% Convert zero rates in percentage
zRates = zRates * 100;

end
