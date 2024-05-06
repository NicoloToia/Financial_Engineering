function X = Trinomial_Tree(dates,a,sigma)

    % time step
    dt = 1;

    % calculate d_r
    d_r = sigma*sqrt(3*(1-exp(-2*a*dt))/(2*a))

    % calculate u_hat
    mu_hat = 1- exp(-a*dt)

    % calculate l_max
    l_max = ceil((1-sqrt(2/3))/mu_hat)
    l_min = -l_max

    % find the levels of the trinomial tree
    l= [0:1:l_max];


    % compute the probabilities of the trinomial tree

    % compute scheme A probabilities
    pu_A = 1/2*(1/3 - mu_hat.*l + mu_hat^2.*l.^2);
    pd_A = 1/2*(1/3 + mu_hat.*l + mu_hat^2.*l.^2);
    ps_A = 2/3 - mu_hat^2.*l.^2; % prob of staying at the same level

    % compute scheme C probabilities
    ps_C = 1/2*(7/3 - 3*mu_hat*l_max + mu_hat^2*l_max^2);
    pdd_C = 1/2*(1/3 - mu_hat*l_max + mu_hat^2*l_max^2);
    pd_C = -1/3 + 2*mu_hat*l_max - mu_hat^2*l_max^2;

    % compute scheme B probabilities
    pu_B = 1/2*(1/3 - mu_hat*l_min + mu_hat^2*l_min^2);
    pd_B = 1/2*(1/3 + mu_hat*l_min + mu_hat^2*l_min^2);
    ps_B = 2/3 - mu_hat^2*l_min^2;
     
    N = length(dates);

    leaves_Tree = 0*ones(1,2*l_max);

    % for i = N-1:-1:1
        

    % end
    X = 1;
end