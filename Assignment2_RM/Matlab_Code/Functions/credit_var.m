
function [VaR_DO, VaR_DM] = credit_var(N, N_issuers, Q, FV, E_FV, rho)
    %% Barriers and simulated P/L for a single IG issuer
    bD = norminv(Q(1,3));                       % Barrier to Def.
    bd = norminv(Q(1,2)+Q(1,3));                % Barrier to down (HY)
    bu = 1.0e6;                                     % Barrier to up (never)
    L_D = (FV(3) - E_FV)/N_issuers;                 % Loss if default
    L_d = (FV(2) - E_FV)/N_issuers;                 % Loss if down
    L_i = (FV(1) - E_FV)/N_issuers;                 % Profit if unchanged
    L_u = 0;                                        % Profit if up 

    %% Monte Carlo
    % Generate N standard 1-d normal variables (macro-factor)
    Y = randn(N,1);                      % Realization of the single factor

    %% Idiosynchratic Monte Carlo
    % zeros(N,1);        % double: Realization of the AVR 
    scen_D = zeros(N,1);   % int.: Count of default events in each MC scenario
    scen_d = zeros(N,1);   % int.: Count of "down" events in each MC scenario
    scen_i = zeros(N,1);   % int.: Count of "identical rating" events in each MC scenario
    scen_u = zeros(N,1);   % int.: Count of "up" events in each MC scenario
    h = waitbar(0, 'Loading...'); % Initialize waitbar
    for i = 1:N_issuers
        V = rho * Y + sqrt(1 - rho^2) * randn(N,1); % Realization of the idiosyncratic factor
        % Update the event counts
        scen_D = scen_D + (V < bD);
        scen_d = scen_d + (V < bd & V >= bD);
        scen_i = scen_i + (V >= bd & V < bu);
        scen_u = scen_u + (V >= bu);
        % Update waitbar
        waitbar(i/N_issuers, h, sprintf('Loading... %d/%d', i, N_issuers));
    end

    close(h); % Close waitbar after loop completes
    % test on marginal distributions
    D_count = sum(scen_D)/N;        % Average number of default events in portfolio
    d_count = sum(scen_d)/N;        % Average number of "down" events in portfolio
    i_count = sum(scen_i)/N;        % Average number of "identical rating" events in portfolio
    u_count = sum(scen_u)/N;        % Average number of "up" events in portfolio
    % Are MC results in agreement with the input rating transition matrix?
    pD_simulated =D_count / N_issuers;
    pd_simulated =d_count / N_issuers;
    pi_simulated =i_count / N_issuers;
    pu_simulated =u_count / N_issuers;

    %% Questions 2 and 3 (and discussion): Credit VaR

    % Default only case (DO)
    EL = L_D * pD_simulated + L_d * pd_simulated + L_i * pi_simulated + L_u * pu_simulated;
    L_y = quantile(-L_D * scen_D,0.999);
    VaR_DO = L_y - EL;

    % Default and migration case (DM)
    L_y = quantile(-L_D * scen_D - L_d * scen_d - L_i * scen_i,0.999);
    VaR_DM = L_y - EL;

end