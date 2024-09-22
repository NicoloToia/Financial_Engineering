function B_contvalue = B_contvalue(l, mu_hat, sigma, sigma_star, a, dt, fwdDF_present, next_value, D_x, x)

% l of the current node
% fwdDF_present right column

B_scheme = [2, 1, 0];

stoch_disc = zeros(3,1);

% do the equation above in a cycle the stoch_disc is a vector 3x1
for i = 1:3
    stoch_disc(i) = fwdDF_present(i) * exp(-0.5 *sigma_star^2 - (sigma_star/sigma) ...
            * (exp(-a*dt) * D_x * B_scheme(i) + mu_hat* x(i)) );
end

[p_u, p_m, p_d] = schemeB(l, mu_hat);
probabilities = [p_u, p_m, p_d]';

B_contvalue =  sum(next_value(end-2:end) .* probabilities .* stoch_disc);

end