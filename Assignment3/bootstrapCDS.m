function [datesCDS, survProbs, intensities] =  bootstrapCDS(datesDF, discounts, datesCDS, spreadsCDS, flag, recovery)
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
EU_ACT_365 = 3;
EU_30_360 = 6;
term1=0;
term2=0;
if flag == 1
    % settlment
    survProbs(1)=1;
    % first year
    survProbs(2)=(1-recovery)*survProbs(1)/(spreadsCDS(1)*yearfrac(datesCDS(1),datesCDS(2),EU_30_360)+(1-recovery));

    deltas = yearfrac(datesCDS(1:end-1), datesCDS(2:end), EU_ACT_365);

    for i = 3:length(datesCDS)

        term1 = term1 + (1-recovery)*intExtDF(discounts,datesDF,datesCDS(i-1))*(survProbs(i-2)- survProbs(i-1)) ;

        term2 = term2 + spreadsCDS(i-1)*yearfrac(datesCDS(i-2),datesCDS(i-1),EU_30_360)*intExtDF(discounts,datesDF,datesCDS(i-1))*survProbs(i-1) ;

        term3 = (1-recovery)*intExtDF(discounts,datesDF,datesCDS(i))*survProbs(i-1) ;

        term4 = (spreadsCDS(i-1)*yearfrac(datesCDS(i-2),datesCDS(i-1),EU_30_360)+(1-recovery))*intExtDF(discounts,datesDF,datesCDS(i));

        survProbs(i) = ( term3 + term1 - term2 ) / term4 ;

    end
    intensities= -1 ./ deltas .* log(survProbs(2:end)/survProbs(1:end-1));
    survProbs = survProbs';
end

end
