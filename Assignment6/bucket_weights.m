function [weights_buck, weights_upto50] = bucket_weights(buckets, selected_bucket)
    % bucket_weights: Compute the bucket weights
    %   buckets: The bucket sizes
    %   selected_bucket: The selected bucket
    %   weights: The bucket weights

    %OUTPUT
    %   weights_buck: The bucket weights
    %   weights_upto50: The weights upto 50

    %if selected_bucket is not in the buckets array, return an error
    if ~ismember(selected_bucket, buckets)
        error('The selected bucket is not in the buckets array');
    end
    %find the position of selected bucket in buckets
    pos = buckets==selected_bucket;
    weights_upto50 = zeros(1, buckets(end));
    if pos>1
        preceding_bucket = buckets(pos-1); %find the preceding bucket
        %linear interpolation to find the weights upto 50 between the preceding bucket and the selected bucket
        weights_upto50(preceding_bucket:selected_bucket) = interp1([preceding_bucket, selected_bucket], [0, 1], preceding_bucket:selected_bucket);
    else 
        weights_upto50(1:selected_bucket) = 1;
    end 
    if pos<buckets(end)
        following_bucket = buckets(pos+1); %find the following bucket
        %linear interpolation to find the weights upto 50 between the selected bucket and the following bucket
        weights_upto50(selected_bucket:following_bucket) = interp1([selected_bucket, following_bucket], [1, 0], selected_bucket:following_bucket);
    else
        weights_upto50(selected_bucket:end) = 1;
    end
end