function Phi = Poisson1D(f, Lx)

%---------------------------------------------
%   1-D Fast Fourier periodic Poisson solver
%---------------------------------------------

%
%	Description
%
%	Using 1-D fast Fourier transform to solve
%	1-D poisson equation
%

%
%   Parameters
%
%   f -> Right function column vector
%   Lx -> Length of x side
%

%
%   Acceptable boundary conditions
%
%   Periodic boundary condition
%   (#Right function must be periodic)
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.16
%

N = length(f);

kx = (2*pi/Lx)*[0:(N/2-1) (-N/2):(-1)]; %->Wavenumber space grid

delsq = -kx.^2;	%->Laplacian matrix acting on the wavenumbers
delsq(1) = 1;

%Spectral inversion of Laplacian
f_fft = fft(f);
Phi = real(ifft(f_fft./delsq));
Phi = Phi - Phi(1);
