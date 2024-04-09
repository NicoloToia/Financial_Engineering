function I = FFT(f, M, xi_1)
% I = FFT(integrand, M, xi_1)
%
% This function computes the Fourier Transform of the input integrand
% using the FFT algorithm.
%
% INPUTS:
%   f: function handle to the integrand
%   M: N = 2^M number of points to use in the FFT
%   xi_1: the fourier inferior limit
%
% OUTPUTS:
%   I: The integral of the integrand
%

% compute N
N = 2^M;

% compute the dxi value
d_xi = -2*xi_1/N;
xi = xi_1:d_xi:-xi_1-d_xi;

% compute the x values
dz = 2 * pi / (N * d_xi);
z = -pi/d_xi:dz:pi/d_xi-dz;

% evaluate the integrand at the xi values
f_j = f(xi);

f_tilde_j = f_j .* exp(-1i * (0:N-1) .* d_xi .* z(1));

% compute the omegas
omega = exp(-1i * 2 * pi / N) .^ (0:N-1)' * (0:N-1);

% compute the FFT vector
FFT = omega * f_tilde_j';

% compute the integral by multiplying by the prefactor
I = d_xi * exp(-1i * xi_1 .* z') .* FFT;

end