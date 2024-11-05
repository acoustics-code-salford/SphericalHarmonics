function [R, Q] = GetRotationMatrix(a,b,c,Order)
%
% [R, Q] = GetRotationMatrix(a,b,c,Order)
%
% Computes a rotation matrix R between the coefficients of Spherical
% Harmonic sound field descriptions before and after rotation. Note that R
% is block diagonal, since rotations only involve coefficients within an
% order, so it's returned as a sparse matrix. R is square & size (Order+1)^2.
% 
% Essentially this is equations 3.3.37 and 3.3.39 of:
% Gumerov, N., & Duraiswami, R. (2005). Fast Multipole Methods for the
% Helmholtz Equation in Three Dimensions (1st ed.). Elsevier Science.
%
% The rotation is actually comprised of three rotations, as detailed on
% page 121 and Eq. 3.3.12:
% 1) Rotation by a radians around z axis of the original coordinate system*
% 2) Rotation by b radians around the y axis of the transitionary coordinate system
% 3) Rotation by c radians around the z axis of the new coordinate system
%
% * note that the formulation there actually rotates by pi-a; this script
% makes that substitution so that a = 0 means no rotation (rather more intuitive!)
% 
% Optionally also returns a 3 x 3 matrix Q, which is the tensor product
% between the original and transformed coordinate system unit vectors.
%
% Arguments:
% a     First rotation angle in radians (real scalar)
% b     Second rotation angle in radians (real scalar)
% c     Third rotation angle in radians (real scalar)
% Order Spherical Harmonic Order (non-negative real integer scalar)
%
% This file is part of the supporting data produced to accompany the
% journal article "A framework for auralization of boundary element method
% simulations including source and receiver directivity” by Jonathan A.
% Hargreaves, Luke R. Rendell and Yiu W. Lam, which was submitted for
% publication in the Journal of the Acoustical Society of America on 13th
% August 2018.
%
% The work was supported by the UK Engineering and Physical Sciences
% Research Council [grant numbers EP/J022071/1 and EP/K000012/1 “Enhanced
% Acoustic Modelling for Auralisation using Hybrid Boundary Integral
% Methods”]. It is provided without warrantee under a Creative Commons
% ‘Attribution’ licence (CC BY); see
% http://creativecommons.org/licenses/by/4.0/ for more information.

% Argument checking:
if ~isscalar(a) && any(imag(a(:))~=0)
    error('a must be a real scalar')
end
if ~isscalar(b) && any(imag(b(:))~=0)
    error('b must be a real scalar')
end
if ~isscalar(c) && any(imag(c(:))~=0)
    error('c must be a real scalar')
end
if ~isscalar(Order) || any(Order(:)<0) || any(imag(Order(:))~=0) || any(mod(Order(:),1)~=0)
    error('Order must be a non-negative real integer scalar')
end

% Allocate R:
R = spalloc((Order+1)^2, (Order+1)^2, (Order+1) .* (Order+2) .* (2*Order+3) ./ 6);

% Loop over SH order:
for n = 0:Order
    for m1 = -n:n
        for m2 = -n:n
            
            % Evalute Eq. 3.3.39:
            if m1 > 0
                ep1 = (-1)^m1;
            else
                ep1 = 1;
            end
            if m2 > 0
                ep2 = (-1)^m2;
            else
                ep2 = 1;
            end
            H = 0;
            for s = max(0, -(m1+m2)):min(n-m1,n-m2)                
                H = H + (-1)^(n-s) * cos(b/2)^(2*s+m2+m1) * sin(b/2)^(2*n-2*s-m2-m1)...
                      / (factorial(s) * factorial(n-m1-s) * factorial(n-m2-s) * factorial(m1+m2+s));
            end
            H = H * ep1 * ep2 * sqrt(factorial(n+m2)*factorial(n-m2)*factorial(n+m1)*factorial(n-m1));
            
            % Evaluate Eq. 3.3.37:
            R(sub2indSH(m2,n), sub2indSH(m1,n)) = (-1)^m1 * exp(-1i*m1*a) * exp(-1i*m2*c) * H;
    
        end
    end
end

% Compute Q if required, using Eq. 3.3.12:
if nargout > 1
    Q1 = [sin(a), cos(a), 0; -cos(a), sin(a), 0; 0, 0, 1];
    Q2 = [-1, 0 0; 0, -cos(b), sin(b); 0, sin(b), cos(b)];
    Q3 = [sin(c), cos(c), 0; -cos(c), sin(c), 0; 0, 0, 1];
    Q = Q3 * Q2 * Q1;
end
