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
%   intensities: intensitiesz


term1=0;
term2=0;

if flag == 1
  
    t0 = 733457;

    survProbs(1)=1;
    survProbs(2)=(1-recovery)*survProbs(1)/(spreadsCDS(1)*yearfrac(t0,datesCDS(2))+(1-recovery));

    for i = 3:length(datesCDS)

        term1 = term1 + (1-recovery)*intExtDF(discounts,datesDF,datesCDS(i-1))*(survProbs(i-2)- survProbs(i-1)) ;

        term2 = term2 + spreadsCDS(i-1)*yearfrac(t0,datesCDS(i-1))*intExtDF(discounts,datesDF,datesCDS(i-1))*survProbs(i-1) ;

        term3 = (1-recovery)*intExtDF(discounts,datesDF,datesCDS(i))*survProbs(i-1) ;

        term4 = (spreadsCDS(i-1)*yearfrac(t0,datesCDS(i-1))+(1-recovery))*discounts(i) ;

        survProbs(i) = ( term3 + term1 - term2) / term4 ;

        intensities(i) = -log(survProbs(i))*(datesCDS(i)-t0);
    end
end

end