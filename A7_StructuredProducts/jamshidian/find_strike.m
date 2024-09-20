function strikes = find_strike(alfa,omega,strike, Talfa, node_dates, fwd_discount_nodes, a, sigma, t0)

ACT_365 = 3;

j=1;
for i = alfa+1 :omega

    deltas(j) = yearfrac(Talfa,node_dates(i),ACT_365);
    
    sigma_hjm(j) = sigmaHJM(a, sigma, node_dates(i), Talfa);

    %integral_ZC(j) = IntHJM(a, sigma, Talfa, t0, node_dates(i)-Talfa)

    j = j+1;
end

integral_ZC = integral_jam(t0,Talfa,node_dates(alfa+1:omega),sigma,a)';

B = @(x)  fwd_discount_nodes'.*exp(-x.*sigma_hjm/sigma - 0.5.*integral_ZC);

 coupon = strike*deltas;

 coupon(end) = coupon(end) + 1;

P = @(x) sum(coupon.*B(x));

x_star = fzero( @(x)P(x) - 1, 0);

strikes = B(x_star);


end