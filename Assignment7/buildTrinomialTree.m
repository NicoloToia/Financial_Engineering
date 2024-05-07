function Tree = buildTrinomialTree(a, sigma, dt, ttm)

    % number of time steps
    N = ttm/dt;

    % calculate mean and variance
    mu_hat = 1- exp(-a*dt)
    sigma_hat = sigma * sqrt((1-exp(-2*a*dt))/(2*a))

    % calculate d_r
    d_r = sigma_hat * sqrt(3);


    % calculate l_max and l_min of the tree
    l_max = ceil((1-sqrt(2/3))/mu_hat);

    % initialize the tree
    Tree = cell(2*l_max+1, N);
    r_0 = 0;
    Tree{1} = r_0;

    % create the tree for the interest rate
    for i = 1:N
        % iterate over the elements of the cell
        nodes = Tree{i};
        % create the new nodes
        if length(nodes) < 2 * l_max + 1
            new_nodes = zeros(length(nodes)+2, 1);
            % write the new nodes
            for j = 1:length(nodes)
                new_nodes(j) = nodes(j) + d_r;
                new_nodes(j+1) = nodes(j);
                new_nodes(j+2) = nodes(j) - d_r;
            end
        else
            new_nodes = nodes;
        end

        % update the tree
        Tree{i+1} = new_nodes;

    end

end