function m = Moments1D(f, vx, order)

%-------------------------------------------------
%   Moments of 1-D (1C1V) distribution function
%-------------------------------------------------

%
%	Description
%
%   Calculate the zeroth, first and second order moments
%   of 1-D distribution function by using trapz method
%

%
%   Parameters
%
%   f -> 1-D distribution function
%   vx -> Velocity coordinates
%   order -> Order of moment
%

%
%   Acceptable input function
%
%   A. The second dimension of f must be velocity
%   B. vx must be a row vector
%   C. order is a string, which can take 'zeroth', 'first', 
%      'second' and 'secondm'
%

%
%   Author: Guanshan Pu; Last modified: 2021.04.19
%

switch order
    case 'zeroth'
        
        m = trapz(vx, f, 2);
    
    case 'first'

        vf = bsxfun(@times, vx, f);
        m = trapz(vx, vf, 2);
        m = m./trapz(vx, f, 2);
        
    case 'second'
       
        v2f = bsxfun(@times, vx.^2, f);
        m = trapz(vx, v2f, 2)./trapz(vx, f, 2);
        
    case 'secondm'
        
        v = Moments1D(f, vx, 'first');
        
        [ngx, ~] = size(f);

        vt = repmat(vx, ngx, 1);
        vt = vt - repmat(v, 1, length(vx));

        tf = vt.^2.*f;

        m = trapz(vx, tf, 2);
        m = m./trapz(vx, f, 2);

    otherwise
        error('This program can only calculate 0~2nd moment')
end