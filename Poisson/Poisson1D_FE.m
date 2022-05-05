function Phi = Poisson1D_FE(f, Lx, bc)

%--------------------------------------
%   1-D Finite element Poisson solver
%--------------------------------------

%
%	Description
%
%	Using Galerkin method with 1-D linear 
%	basis to solve 1-D non-periodic Poisson
%	equation
%

%
%   Parameters
%
%   f -> Right function column vector
%   Lx -> Length of system
%   bc -> Boundary conditions
%

%
%	bc is a MATLAB structure and must contains
%	two of following fields:
%
%	(n1, d2) or (d1, d2) or (d1, n2)
%

%
%   Acceptable boundary conditions
%
%   A. left:Dirichlet, right:Dirichlet
%   B. left:Dirichlet, right:Neumann
%   C. left:Neumann, right:Dirichlet
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.16
%

N = length(f);dx = Lx/(N - 1); 

u = ones(N, 1);

%Generate coefficient matrix for linear test function
A = spdiags([-u, 2*u, -u], [-1 0 1], N, N);     %->Laplacian
A = full(A);

B = spdiags([u/6, 2*u/3, u/6], [-1 0 1], N, N); %->Right function
B = full(B);

%Select the boundary conditions
if isfield(bc, 'd1') && isfield(bc, 'd2')     %->left:Dirichlet, right:Dirichlet
    
    %->Value of test function at boundary is zero
    
    A = A(2:end-1, 2:end-1)/dx^2;
    B(1, :) = [];B(end, :) = [];
    
    f_p = B*f;

    K = zeros(N-2, 1);
    K(1) = bc.d1;K(end) = bc.d2;
    K = K/dx^2;
    
    f_p = f_p + K;

    Phi = [bc.d1; A^-1*f_p; bc.d2];

    
elseif isfield(bc, 'd1') && isfield(bc, 'n2') %->left:Dirichlet, right:Neumann
    
    %->Value of test function at left boundary is zero
    
    A(end) = 1;
    A = A(2:end, 2:end)/dx^2;
    B(1, :) = [];B(end) = 1/3;
    
    f_p = B*f;
    
    K = zeros(N-1, 1);
    K(1) = bc.d1/dx^2;K(end) = bc.n2/dx;
    
    f_p = f_p + K;
    
    Phi = [bc.d1; A^-1*f_p];
    
elseif isfield(bc, 'n1') && isfield(bc, 'd2') %->left:Neumann, right:Neumann
    
    %->Value of test function at right boundary is zero
    
    A(1) = 1;
    A = A(1:end-1, 1:end-1)/dx^2;
    B(end, :) = [];B(1) = 1/3;
    
    f_p = B*f;
    
    K = zeros(N-1, 1);
    K(1) = bc.n1/dx;K(end) = bc.d2/dx^2;
    
    f_p = f_p + K;
    
    Phi = [A^-1*f_p; bc.d2];
    
end


