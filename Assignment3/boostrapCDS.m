[datesCDS, survProbs, intensities] =  bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)

%    INPUT:   
%   datesDF: dates of the discount factors
%   discounts: discount factors
%   datesCDS: dates of the CDS
%   spreadsCDS: spreads of the CDS
%   flag: 1 if the CDS spreads are in basis points, 0 if they are in decimals
%   recovery: recovery rate

%   OUTPUT: 
%   datesCDS: dates of the CDS
%   survProbs: survival probabilities
%   intensities: intensities

