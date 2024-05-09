function A_contvalue = A_contvalue(l, mu_hat, sigma, sigma_star, a, dt, fwdDF_present, next_value, D_x, x)

% l of the current node
% fwdDF_present right column

A_scheme = [1, 0, -1];


% % do the equation above in a cycle the stoch_disc is a vector 3x1
% for i = 1:3
%     stoch_disc(i) = fwdDF_present(l-1) * exp(-0.5 *sigma_star^2 - (sigma_star/sigma) ...
%             * (exp(-a*dt) * D_x * A_scheme(i) + mu_hat* x(i)) );
% end

stoch_disc(1) = fwdDF_present(l-1) * exp(-0.5 *sigma_star^2 - (sigma_star/sigma) ...
            * (exp(-a*dt) * D_x * A_scheme(1) + mu_hat* x(l-1)) );
stoch_disc(2) = fwdDF_present(l) * exp(-0.5 *sigma_star^2 - (sigma_star/sigma) ...
            * (exp(-a*dt) * D_x * A_scheme(2) + mu_hat* x(l)) );
stoch_disc(3) = fwdDF_present(l+1) * exp(-0.5 *sigma_star^2 - (sigma_star/sigma) ...
            * (exp(-a*dt) * D_x * A_scheme(3) + mu_hat* x(l+1)) );


[p_u, p_m, p_d] = schemeA(l, mu_hat);
probabilities = [p_u, p_m, p_d]';

A_contvalue =  sum(next_value(1:3) .* probabilities .* stoch_disc');

end