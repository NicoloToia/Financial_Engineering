function FT = MC_NIG(F_0, kappa, t, alpha, eta, sigma, x, N, discount_1y)

% compute the Laplace exponent
ln_L = @(omega) t/kappa * (1 - alpha)/alpha * ...
    (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );
% draw the standard normal random variables
g = randn(N, 1);
% draw the inverse gaussian random variables
G = random('inversegaussian', 1, t/kappa, N, 1);

ft = sqrt(t) * sigma * sqrt(G) .* g - (0.5 + eta) * t * sigma^2 * G - ln_L(eta);

FT = F_0 * exp(ft);

end