function bucket_sensitivities = coarseBuckets(sens_dates, sensitivities, bucket_end_dates)

% initialize the sensitivities
bucket_sensitivities = zeros(4, 1);

bucket_sensitivities(1) = bucketWeights(sensitivities, sens_dates, 0, 2, 5, true, false);

bucket_sensitivities(2) = bucketWeights(sensitivities, sens_dates, 2, 5, 10, false, false);

bucket_sensitivities(3) = bucketWeights(sensitivities, sens_dates, 5, 10, 15, false, false);

bucket_sensitivities(4) = bucketWeights(sensitivities, sens_dates, 10, 15, 50, false, true);

end