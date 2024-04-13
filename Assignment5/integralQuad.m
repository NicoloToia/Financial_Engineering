function I = integralQuad(phi, M, queryPoints)
% Compute the price of the integral using the quadrature method between lower and upper
%
% Inputs:
%   phi: function to integrate
%   M: N = 2^M, number of nodes for the quadrature
%   queryPoints: points to evaluate the integral
%   lower: lower bound of the integral
%   upper: upper bound of the integral
%
% Outputs:
%   I: integral of the function phi between lower and upper
%

% Transform the function to match the Lewis definition
f = @(xi,x) exp(-1i * xi .* x) / (2 * pi) .* phi(-xi-1i/2) .* 1 ./ (xi.^2 + 1/4);

% find the extrema of the integral if not provided
lower = 0;
while real(f(lower, 0)) > 1e-10
    lower = lower - 1;
end

upper = 0;
while real(f(upper, 0)) > 1e-10
    upper = upper + 1;
end

% Compute the nodes of integration
xi = linspace(lower, upper, 2^M);

% Compute the weights
w = ones(1, 2^M);
w(1) = 1/2;
w(end) = 1/2;

% Compute the integrals for each query point
I = zeros(size(queryPoints));
for i = 1:length(queryPoints)
    I(i) = sum(w .* f(xi, queryPoints(i)));
end

% return only the real part
I = real(I);

end