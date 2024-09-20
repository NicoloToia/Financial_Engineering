function fwd_spot_node = compute_fwdSpot(dates, node_dates, discounts, N_step)
% compute_fwdSpot: computes the forward discount for each node of the tree

% compute the discount factors at the nodes
discount_nodes = zeros(N_step +1, 1);
discount_nodes(1) = 1;

discount_nodes(2:end) = intExtDF(discounts, dates, node_dates(2:end));


% compute the forward discount in the nodes 
fwd_spot_node =  discount_nodes(2:end) ./ discount_nodes(1:end-1);

end