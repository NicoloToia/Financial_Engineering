function prob = closed_formula_digital(S,K,T,r,q,eta,alpha,sigma,kappa)

    k= log(S/K) + (r - q)*T;

    ln_L = @(omega) T/kappa * (1 - alpha)/alpha * ...
        (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );

    ln_L_eta = ln_L(eta);

    phi = @(xi) exp(-1i * xi * ln_L_eta) .* exp( ln_L (0.5 * ((xi.^2) + 1i * (1+2*eta) .* xi)));

    integrand = @(u) real(exp(1i*u*k).*phi(u)./(1i*u));

    I = integral(integrand, 0, Inf);

    prob = 0.5 - 1/pi * I;

end