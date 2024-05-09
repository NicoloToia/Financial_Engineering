function FT = MC_NIG(F_0, alpha, sigma, kappa, eta, t,  N)
% 
%INPUTS: 
%
%   F_0: initial price of the underlying corresponding to the forward
%   alpha = parameter to identify the method (in this case, alpha = 0.5)
%   sigma: volatility of the underlying
%   kappa: vol of vol
%   eta: skewness
%   t = yearfrac considered 
%   N = number of simulations for MonteCarlo 

%OUTPUT: 
% 
%   FT = dynamic of the underlying of the forward simulated through MC

%compute the Laplace Exponent
ln_L = @(omega) t/kappa * (1 - alpha)/alpha * ...
    (1 - (1 + (omega .* kappa * sigma^2)/(1-alpha)).^alpha );
% draw the standard normal random variables
g = randn(N, 1);
% draw the inverse gaussian random variables
G = random('inversegaussian', 1, t/kappa, N, 1);

ft = sqrt(t) * sigma * sqrt(G) .* g - (0.5 + eta) * t * sigma^2 * G - ln_L(eta);
if length(F_0) == 1 
    FT = F_0 .* exp(ft);
else 
    %in this case we have already simulated the first part of the process,
    %we want to simulate what happens between first and second year
    FT = zeros(N, 1); 
    for ii = 1:N 
        FT(ii, 1) = F_0(ii, 1) * exp(ft(ii));
    end 
end 
