function bucket_sensitivities = deltaCoarseBuckets(settlement_date, sens_dates, sensitivities)
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
bucket_sensitivities = zeros(4, 1);

% compute the bucket sensitivities for the first bucket
bucket_sensitivities(1) = bucketWeights(settlement_date, sensitivities, sens_dates, 0, 2, 5, true, false);

% compute the bucket sensitivities for the second bucket
bucket_sensitivities(2) = bucketWeights(settlement_date, sensitivities, sens_dates, 2, 5, 10, false, false);

% compute the bucket sensitivities  for the third bucket
bucket_sensitivities(3) = bucketWeights(settlement_date, sensitivities, sens_dates, 5, 10, 15, false, false);

% compute the bucket sensitivities  for the fourth bucket
bucket_sensitivities(4) = bucketWeights(settlement_date, sensitivities, sens_dates, 10, 15, 50, false, true);

end