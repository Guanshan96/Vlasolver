function grad = Gradient1D(f, Lx)

%---------------------------------------------
%   1-D Fast Fourier numerical gradient
%---------------------------------------------

%
%   Parameters
%
%   f -> Input function
%   Lx -> Length of x side
%

%
%   Acceptable input function
%
%   Periodic function
%	(#f must be periodic)
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.16
%

N = length(f);

kx = (2*pi/Lx)*[0:(N/2-1) (-N/2):(-1)]; %Wavenumber vector

f_fft = fft(f);
grad = real(ifft(1i*kx.*f_fft));