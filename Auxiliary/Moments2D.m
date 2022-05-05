function m = Moments2D(f, vx, vy, order)

[ngy, ngx, ~, ~] = size(f);

switch order
    case 'zeroth'
        
        m = trapz(vx, f, 4);
        m = trapz(vy, m, 3);
    
    case 'first'
        
        m = zeros(ngy, ngx, 2);
        
        vf = trapz(vy, f, 3);
        vf = bsxfun(@times, reshape(vx, [1 1 l length(vx)]), vf);
        
        m(:, :, 1) = trapz(vx, vf, 4);
        
        vf = trapz(vx, f, 4);
        vf = bsxfun(@times, reshape(vy, [1 1 length(vy) 1]), vf);
        
        m(:, :, 2) = trapz(vy, vf, 3);
        
    case 'second'
       
        m = zeros(ngy, ngx, 3);
        
        v2f = trapz(vy, f, 3);
        v2f = bsxfun(@times, reshape(vx.^2, [1 1 1 length(vx)]), v2f);
        
        m(:, :, 1) = trapz(vx, v2f, 4);
        
        v2f = trapz(vx, f, 4);
        v2f = bsxfun(@times, reshape(vy.^2, [1 1 length(vy) 1]), v2f);
        
        m(:, :, 2) = trapz(vy, v2f, 3);
        
        v2f = bsxfun(@times, reshape(vx, [1 1 1 length(vx)]), f);
        v2f = bsxfun(@times, reshape(vy, [1 1 length(vy) 1]), v2f);
        
        v2f = trapz(vx, v2f, 4);
        m(:, :, 3) = trapz(vy, v2f, 3);
        
    case 'secondm'
        
        v = Moments2D(f, vx, vy, 'first');
        
        [ngx, ~] = size(f);

        vt = repmat(vx, ngx, 1);
        vt = vt - repmat(v, 1, length(vx));

        tf = vt.^2.*f;

        m = trapz(vx, tf, 2);
        m = m./trapz(vx, f, 2);

    otherwise
        error('')
end