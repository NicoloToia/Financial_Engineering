function price = bermudanSwaption(a, sigma, dt, ttm, reset_dates, discounts, dates)

% build the trinomial tree
[l_max, mu, trinomial_tree] = buildTrinomialTree(a, sigma, dt, ttm);

% compute the discount factors at the reset dates
reset_DF = intExtDF(discounts, dates, reset_dates);

% continuation values
values = zeros(size(Tree{length(Tree)}));

% iterate over the tree
for i = length(Tree)-1:-1:l_max+1

    % iterate over the nodes
    nodes = Tree{i};
    continuation_values = zeros(size(nodes));

    for j = 1:nodes

        % compute the forward discount factor
        B_dt = 

        % apply the scheme and compute the intrinsic values
        if j == 1 % scheme C
            % compute the probabilities of the trinomial tree
            [p_u, p_m, p_d] = schemeC(l_max, mu);
            % compute the stochastic discount factors
            [D_u, D_m, D_d]
            % compute the mean of the first three elements of the past nodes
            continuation_values(j) = p_u * values(j) + p_m * values(j+1) + p_d * values(j+2);
        elseif j == length(nodes) % scheme B
            [p_u, p_m, p_d] = schemeB(-l_max, mu);
            continuation_values(j) = p_u * values(j-2) + p_m * values(j-1) + p_d * values(j);
        else % scheme A
            [p_u, p_m, p_d] = schemeA(l_max - j + 1, mu);
            continuation_values(j) = p_u * values(j-1) + p_m * values(j) + p_d * values(j+1);
        end

        % compute the intrinsic value
        ttm_swap = dt * (length(Tree) - i);
        ttm_option = dt * i;
        % compute the swap rate
        swap_rate = swapRateHW(nodes(j), fwd_DF(i), sigma, a, tau);

    end
    
end
