function Phi = Poisson1D_FD(f, Lx, bc)

%--------------------------------------
%   1-D Finite difference Poisson solver
%--------------------------------------

%
%	Description
%
%	Using second order finite difference method
%	to solve 1-D non-periodic Poisson equation
%

%
%   Parameters
%
%   f -> Right function column vector
%   Lx -> Length of system
%   bc -> Boundary conditions
%

%
%	bc is a string which value is 'periodic'
%	when the boundary condition is periodic
%
%	bc is a MATLAB structure and must contains
%	two of following fields:
%
%	(n1, d2) or (d1, d2) or (d1, n2)
%
%	When non-periodic is applied
%

%
%   Acceptable boundary conditions
%
%   A. left:Dirichlet, right:Dirichlet
%   B. left:Dirichlet, right:Neumann
%   C. left:Neumann, right:Dirichlet
%   D. Periodic
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.18
%

N = length(f); dx = Lx/(N - 1);
un = ones(N, 1);

A = spdiags([-un, 2*un, -un], [-1 0 1], N, N);
A = full(A);

if strcmp(bc, 'periodic') %->Periodic
    
    f_p = f(1:end - 1);

    A = A(2:end, 2:end);
    A(1, end) = -1;A(end, 1) = -1;A(end, :) = 1;
    A = A^-1;

    f_p(end) = 0;

    Phi = A*(f_p*dx^2);Phi = [Phi; Phi(1)];
    
elseif isfield(bc, 'd1') && isfield(bc, 'd2') %->left:Dirichlet, right:Dirichlet
    
    A = A(2:end-1, 2:end-1);
    A = A^-1;

    f_p = f(2:end-1);
    f_p(1) = f_p(1) + bc.d1/dx^2;
    f_p(end) = f_p(end) + bc.d2/dx^2;

    Phi = [bc.d1; A*(f_p*dx^2); bc.d2];
    
elseif isfield(bc, 'd1') && isfield(bc, 'n2') %->left:Dirichlet, right:Neumann
        
    A = A(2:end, 2:end);
    A(end) = 1;
    A = A^-1;

    f_p = f(2:end);
    f_p(1) = f_p(1) + bc.d1/dx^2;
    f_p(end) = bc.n2/dx;

    Phi = [bc.d1; A*(f_p*dx^2)];
    
elseif isfield(bc, 'n1') && isfield(bc, 'd2') %->left:Neumann, right:Neumann
    
    A = A(2:end, 2:end);
    A(1) = -1;A(1, 2) = 1;
    A = A^-1;

    f_p = f(1:end-1);
    f_p(1) = bc.n1/dx;
    f_p(end) = f_p(end) + bc.d2/dx^2;

    Phi = [A*(f_p*dx^2); bc.d2];
    
end

