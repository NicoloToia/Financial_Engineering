function discount = discount_reset(sigma,a,reset_dates,datenodes,trinomial_tree,fwd_discount_reset,l_max,N_step,dates,ttm)

k = 1; % reset date counter

date_t0 = datenum(dates(1));

discount = zeros(2*l_max+1,N_step);

for i = 1:N_step

    if datenodes(i) == reset_dates(k)
        if i < l_max
            for j = 1 : i
                sigma_hjm = sigma/a * (1 - exp(-a .* (ttm-k)));
                integral_ZC = (sigma/a)^2 *((reset_dates(k) - date_t0) - 3/(2*a) + 1/a * exp(-a* ...
                    (reset_dates(k) - date_t0)) * (2 - 0.5*exp(-a*(reset_dates(k) - date_t0))));
                % up case 
                discount(l_max +1 + j,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max +1+ j, i)*sigma_hjm - 0.5*integral_ZC);

                % down case
                discount(l_max +1- j,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max+1 - j, i)*sigma_hjm - 0.5*integral_ZC);

                % stay case
                discount(l_max+1,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max+1, i)*sigma_hjm - 0.5*integral_ZC);
            end
        else
            for j = 1 : l_max
                sigma_hjm = sigma/a * (1 - exp(-a .* (ttm-k)));
                integral_ZC = (sigma/a)^2 *((reset_dates(k) - date_t0) - 3/(2*a) + 1/a * exp(-a* ...
                    (reset_dates(k) - date_t0)) * (2 - 0.5*exp(-a*(reset_dates(k) - date_t0))));

                    discount(l_max +1 + j,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max +1+ j, i)*sigma_hjm - 0.5*integral_ZC);

                    % down case
                    discount(l_max +1 - j,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max+1 - j, i)*sigma_hjm - 0.5*integral_ZC);
    
                    % stay case
                    discount(l_max + 1,i) = fwd_discount_reset(k) * exp(-trinomial_tree(l_max+1, i)*sigma_hjm - 0.5*integral_ZC);
            end
        end
        k=k+1;
    end
end

end