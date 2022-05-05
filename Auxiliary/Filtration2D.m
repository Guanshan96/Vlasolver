function fk = Filtration2D(fd, Vm, order, bw)

%-------------------------------------
%   2-D filamentation filtration
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

[~, ~, Nvy, Nvx] = size(fd);

%--->Generating filter
if mod(Nvy, 2) == 0
    kvx = (2*pi/(2*Vm(1)))*[0:(Nvx/2-1) (-Nvx/2):(-1)];
else
    kvx = (2*pi/(2*Vm(1)))*[0:(Nvx/2-1) (-Nvx/2):(0)];
end

if mod(Nvx, 2) == 0
    kvy = (2*pi/(2*Vm(2)))*[0:(Nvy/2-1) (-Nvy/2):(-1)];
else
    kvy = (2*pi/(2*Vm(2)))*[0:(Nvy/2-1) (-Nvy/2):(0)];
end

[Kvx, Kvy] = meshgrid(kvx, kvy);

Kvx = reshape(Kvx, [1 1 Nvy Nvx]);Kvy = reshape(Kvy, [1 1 Nvy Nvx]);

filter = exp(-(Kvx/bw(1)).^(2*order(1))-(Kvy/bw(2)).^(2*order(2)));

%--->Performing 2-D FFT on input distribution
fk = fft(permute(fft(fd, [], 3), [1 2 4 3]), [], 3);
fk = permute(fk, [1 2 4 3]);

%--->Filtering
fk = bsxfun(@times, filter, fk);

%--->Performing invert 2-D FFT on filtered distribution
fk = ifft(permute(ifft(fk, [], 3), [1 2 4 3]), [], 3);
fk = real(permute(fk, [1 2 4 3]));