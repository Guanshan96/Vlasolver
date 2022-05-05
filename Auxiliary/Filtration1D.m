function fk = Filtration1D(fd, Vm, order, bw)

%-------------------------------------
%   1-D filamentation filtration
%-------------------------------------

%
%	Description
%
%   1-D filamentation filtration by fast Fourier transform.
%   The second dimension of fd will be filted  
%

%
%   Parameters
%
%   fd -> 2-variable scalar value function
%	Vm -> Maximum velocity
%	order -> Order of filter kernel
%   bw -> Bandwidth of filter kernel
%

%
%   Acceptable input parameter
%   
%   
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.22
%

[~, Nv] = size(fd);Vrange = Vm(2) - Vm(1);

if mod(Nv, 2) == 0
    kx = (2*pi/Vrange)*[0:(Nv/2-1) (-Nv/2):(-1)];
else
    kx = (2*pi/Vrange)*[0:(Nv/2-1) (-Nv/2):(0)];
end
fk = fft(fd, [], 2);

filter = exp(-(kx/bw).^(2*order));

fk = real(ifft(bsxfun(@times, filter, fk), [], 2));
