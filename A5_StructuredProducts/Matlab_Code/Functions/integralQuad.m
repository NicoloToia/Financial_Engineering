function I = integralQuad(phi, M, dz, queryPoints)
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

% compute N
N = 2^M;

% compute the dxi value
d_xi = 2 * pi / (N * dz);
xi_1 = -(N-1)/2 * d_xi;

% check that -xi_1 has a value less than 1e-10
if (f(xi_1, 0) > 1e-10)
    error('The function is not integrable at -xi_1, increase the number of nodes M');
end

% return only the real part
I = zeros(size(queryPoints));

for i = 1:length(queryPoints)
    I(i) = quadgk(@(xi) real(f(xi, queryPoints(i))), xi_1, -xi_1, 'MaxIntervalCount', N);
end

end