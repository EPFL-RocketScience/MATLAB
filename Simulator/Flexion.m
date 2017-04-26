function [T, M] = Flexion(m, l, xs, x, a, Ltot)
%FLEXION Calculate shear and flexion moment in beam according to given mass
%distribution
%   [T, M] = Flexion(m, l, xs, x, a) Calculate shear T and flexion M given a 
%   certain mass m, distributed on a certain length l, starting from xs and 
%   evaluated in x. a is the lateral acceleration and Ltot the total length
%   of the beam.  
%   All values can be vectors and must have the same length. m, l, xs must
%   be column vectors. x must be a line vector.

    % check vector dimensions
    if(~iscolumn(m))
        error('m must be a column vector.')
    elseif(~iscolumn(l))
        error('l must be a column vector.')
    elseif(~iscolumn(xs))
        error('xs must be a column vector.')
    elseif(iscolumn(x))
        error('x must be a line vector.')
    end
    
    % number of masses
    n = length(m);
    
    % number of requested data points
    x_length = length(x);
    
    % distributed force
    p = m./l*a;
    
    % Compute matrices
    P   =   repmat(p, 1, x_length);
    L   =   repmat(l, 1, x_length);
    XS  =   repmat(xs, 1, x_length);
    X   =   repmat(x, n, 1);
    
    % reaction force at x = 0
    AY = P.*L.*(1-(XS+L/2)/Ltot);
    
    % Initialize traction and Flexion matrices
    T = zeros(n, x_length);
    M = zeros(n, x_length);
    
    % Compute traction and flexion for 0<x<xs
    X1      = X<XS;
    T(X1)   = AY(X1);
    M(X1)   = X(X1).*AY(X1);
    
    % Compute traction and flexion for xs<x<xs+l
    X2      = (XS<=X & X<=(XS+L));
    T(X2)   = AY(X2) - (X(X2)-XS(X2)).*P(X2);
    M(X2)   = X(X2).*AY(X2)-(X(X2)-XS(X2)).^2.*P(X2)/2;
    
    % Compute tration and flexion for xs+l<x<Ltot
    X3      = (XS+L)< X;
    T(X3)   = AY(X3) - L(X3).*P(X3);
    M(X3)   = X(X3).*AY(X3) - (X(X3)-XS(X3)-L(X3)/2).*L(X3).*P(X3);
    
    % sum contributions
    T = sum(T, 1);
    M = sum(M, 1);
end

