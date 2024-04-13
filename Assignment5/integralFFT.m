function I = integralFFT(phi, M, log_moneyness, xi_1) 

% Compute the number of nodes and the z for the FFT
N = 2^M;
z_1 = -(N-1)/2*dz;
z = z_1:dz: -z_1; 

% Compute the nodes of integration
dx = 2*pi/(N*dz);
x_1 = -(N-1)/2*dx;
x = x_1:dx:-x_1;


% Compute the mu coefficient for the Martingale condition
mu = -log(L(eta));

% Perform FFT
pref = dx * exp(-1i*x(1)*z); % Computation of the prefactor

% Computation of the characteristic function of the forward
phi = @(x) exp(1i*x*mu) .* L( 1i*x*(1/2+eta) +0.5*x.^2);

% Adjust the f to exploit FFT
f = phi(-x-1i/2)./(x.^2+0.25)/(2*pi);
f = f .* exp(-1i*z_1*dx*(0:N-1));

% Exploit FFT, multiply by the prefactor and the exponential term 
% of moneyness as in Lewis formula
integral = exp(-z/2).*fft(f).*pref ;

% I want only the values that I require and take the real part only
integral = real(interp1(z,integral, moneyness));
prices = B*F0*(1 - integral);

end