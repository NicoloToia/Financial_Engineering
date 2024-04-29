function bucket_sensitivities = vegaCoarseBuckets(settlement_date, sens_dates, sensitivities)
% DELTACOARSEBUCKETS computes the coarse-grained delta-bucket sensitivities for the certificate
%
% INPUTS
%   settlement_date: settlement date
%   sens_dates: dates of the sensitivities
%   sensitivities: sensitivities of the certificate
% OUTPUTS
%   bucket_sensitivities: coarse-grained delta-bucket sensitivities for the certificate
%   (0-2y, 2-5y, 5-10y, 10-15y, 15-50y)

% initialize the sensitivities
bucket_sensitivities = zeros(2, 1);

% compute the bucket sensitivities
bucket_sensitivities(1) = bucketWeights(settlement_date, sensitivities, sens_dates, 0, 5, 15, true, false);

% compute the bucket sensitivities
bucket_sensitivities(2) = bucketWeights(settlement_date, sensitivities, sens_dates, 5, 15, 30, false, true);

end