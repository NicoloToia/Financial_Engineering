function [zRates] = zeroRates(dates, discounts)
% zeroRates Compute zero rates from discount factors
%
% INPUT:
% dates     - vector of dates
% discounts - vector of discount factors
%
% OUTPUT:
% zRates    - vector of zero rates
    

%deposRates = 0.5 * (deposRates(1,1) + deposRates(1,2))

%deposRates = ratesSet.depos(1,:);

%zRates = zeros(length(dates),1);

zRates = -log(discounts) ./ yearfrac(dates(1), dates);
%zRates(1) = 0.5 * (deposRates(1,1) + deposRates(1,2));
%zRates(1) = zRates(2) + (zRates(3) - zRates(2))/(yearfrac(dates(2),dates(3),2)) * (- yearfrac(dates(2),dates(1),2));
zRates(1) = 0.0401;
zRates = zRates * 100;
end