function Si = Spinterp_Natural(f, dx, dt, F, dim)
%-------------------------------------------------------------
%   1-D Parallel cubic spline interpolation with natural BC
%-------------------------------------------------------------

%
%	Description
%
%	A parallel cubic spline interpolation code which
%	support up to 5-D array and use natural boundary
%	condition (Second order derivative at boundary is
%	zero). Only available on NVIDIA GPU
%

%
%   Parameters
%
%   f -> Multi-variable scalar value function
%   dx -> Grid size
%	dt -> Time step
%	F -> Advection speed
%	dim -> Direction of interpolation
%

%
%   Acceptable input parameter
%
%	A. Courant coefficient (Fdt/dx) must smaller than 1
%	B. Size of f along other dimension should equal with 
%	size of F
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.18
%

%----Permute selected dimension to the first dimension----
[n1, n2, n3, n4, n5] = size(f);

switch dim
    case 1
                n = n1;
                d = permute(f, [1 2 3 4 5]);
                d = reshape(d, [n1 n2*n3*n4*n5]);
                
                F = permute(F, [1 2 3 4 5]);
                F = reshape(F, [1 n2*n3*n4*n5]);
    case 2
                n = n2;
                d = permute(f, [2 1 3 4 5]);
                d = reshape(d, [n2 n1*n3*n4*n5]);
        
                F = permute(F, [2 1 3 4 5]);
                F = reshape(F, [1 n1*n3*n4*n5]);
    case 3
                n = n3;
                d = permute(f, [3 2 1 4 5]);
                d = reshape(d, [n3 n2*n1*n4*n5]);
        
                F = permute(F, [3 2 1 4 5]);
                F = reshape(F, [1 n2*n1*n4*n5]);                
    case 4
                n = n4;
                d = permute(f, [4 2 3 1 5]);
                d = reshape(d, [n4 n2*n3*n1*n5]);
        
                F = permute(F, [4 2 3 1 5]);
                F = reshape(F, [1 n2*n3*n1*n5]);
    case 5
                n = n5;
                d = permute(f, [5 2 3 4 1]);
                d = reshape(d, [n5 n2*n3*n4*n1]);
        
                F = permute(F, [5 2 3 4 1]);
                F = reshape(F, [1 n2*n3*n4*n1]);
    otherwise
        error('Dimension cannot be negative or exceeds 5')
end

%----Construct coefficient matrix----
un = ones(n,1)*dx;
A = spdiags([un, 4*un, un], [-1 0 1], n-2, n-2);
A = full(A);A = A^-1;

%----Construct right column vector----
d = (d(3:end, :) - 2*d(2:end-1, :) + d(1:end-2, :))/dx;

%----Calculate interpolation coefficients----
[~, NN] = size(d);
d = A*d*6;d = [zeros(1, NN);d;zeros(1, NN)];

%----Calculate values at query points----
zeta = abs(F)*dt;
zeta_m = dx - abs(F)*dt;

zeta = zeta.^3;
zeta_m = zeta_m.^3;

Si = zeros(size(d));

sig = F > 0;
Si(2:end,sig) = bsxfun(@times, zeta_m(sig), d(2:end, sig))/(6*dx);
Si(2:end,sig) = Si(2:end,sig) - bsxfun(@times, zeta_m(sig), d(1:end-1, sig))/(6*dx);
sig = F <= 0;
Si(1:end-1,sig) = bsxfun(@times, zeta(sig), d(2:end, sig))/(6*dx);
Si(1:end-1,sig) = Si(1:end-1,sig) - bsxfun(@times, zeta(sig), d(1:end-1, sig))/(6*dx);

zeta = zeta.^(2/3);
zeta_m = zeta_m.^(2/3);

sig = F > 0;
Si(2:end,sig) = Si(2:end,sig) + bsxfun(@times, zeta_m(sig), d(1:end-1, sig)/2);
sig = F <= 0;
Si(1:end-1,sig) = Si(1:end-1,sig) + bsxfun(@times, zeta(sig), d(1:end-1, sig)/2);

zeta = zeta.^(1/2);
zeta_m = zeta_m.^(1/2);
sig = F > 0;
Si(2:end,sig) = Si(2:end,sig) - bsxfun(@times, zeta_m(sig), dx*d(1:end-1, sig)/3);
Si(2:end,sig) = Si(2:end,sig) - bsxfun(@times, zeta_m(sig), dx*d(2:end, sig)/6);
sig = F <= 0;
Si(1:end-1,sig) = Si(1:end-1,sig) - bsxfun(@times, zeta(sig), dx*d(1:end-1, sig)/3);
Si(1:end-1,sig) = Si(1:end-1,sig) - bsxfun(@times, zeta(sig), dx*d(2:end, sig)/6);

switch dim
    case 1
                d = permute(f, [1 2 3 4 5]);
                d = reshape(d, [n1 n2*n3*n4*n5]);
    case 2
                d = permute(f, [2 1 3 4 5]);
                d = reshape(d, [n2 n1*n3*n4*n5]);
    case 3
                d = permute(f, [3 2 1 4 5]);
                d = reshape(d, [n3 n2*n1*n4*n5]);
    case 4
                d = permute(f, [4 2 3 1 5]);
                d = reshape(d, [n4 n2*n3*n1*n5]);
    case 5
                d = permute(f, [5 2 3 4 1]);
                d = reshape(d, [n5 n2*n3*n4*n1]);
end

sig = F > 0;
Si(2:end,sig) = Si(2:end,sig) + d(1:end-1, sig);
Si(2:end,sig) = Si(2:end,sig) + bsxfun(@times, zeta_m(sig), d(2:end, sig))/dx;
Si(2:end,sig) = Si(2:end,sig) - bsxfun(@times, zeta_m(sig), d(1:end-1, sig))/dx;
Si(1, sig) = d(1, sig);
sig = F <= 0;
Si(1:end-1,sig) = Si(1:end-1,sig) + d(1:end-1, sig);
Si(1:end-1,sig) = Si(1:end-1,sig) + bsxfun(@times, zeta(sig), d(2:end, sig))/dx;
Si(1:end-1,sig) = Si(1:end-1,sig) - bsxfun(@times, zeta(sig), d(1:end-1, sig))/dx;
Si(end, sig) = d(end, sig);


switch dim
    case 1
        Si = reshape(Si, [n1 n2 n3 n4 n5]);
        Si = permute(Si, [1 2 3 4 5]);
    case 2
        Si = reshape(Si, [n2 n1 n3 n4 n5]);
        Si = permute(Si, [2 1 3 4 5]);
    case 3
        Si = reshape(Si, [n3 n2 n1 n4 n5]);
        Si = permute(Si, [3 2 1 4 5]);
    case 4
        Si = reshape(Si, [n4 n2 n3 n1 n5]);
        Si = permute(Si, [4 2 3 1 5]);
    case 5
        Si = reshape(Si, [n5 n2 n3 n4 n1]);
        Si = permute(Si, [5 2 3 4 1]);
end
