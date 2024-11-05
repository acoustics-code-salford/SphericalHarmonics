function [Y, dY_dbeta, dY_dalpha] = SphericalHarmonic(m,n,alpha,sinbeta,cosbeta)
%
% [Y, dY_dbeta, dY_dalpha] = SphericalHarmonic(m,n,alpha,sinbeta,cosbeta)
% [Y, dY_dbeta, dY_dalpha] = SphericalHarmonic(m,n,alpha,beta)
%
% Computes a Spherical Harmonic function of order (m,n), and it's angular
% derivatives if required.
%
% Arguments:
% r is radius
% alpha is azimuth angle (angle in radians from the positive x axis, with
% rotation around the positive z axis according to the right-hand screw rule)
% beta is polar angle (angle in radians from the positive z axis). It can
% alternatively be specified as two arrays of its cos and sin values. 
% These should all be arrays of identical size.
% m and n should be integer scalars; n should be non-negative and m should be in the range -n<=m<=n
%
% Returned arrays will have numel(alpha) rows and (Order+1)^2 columns.
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
if ~isscalar(n) || ~all(isreal(n(:))) || ~all(mod(n(:),1)==0) || any(n(:)<0)
    error('n must be a non-negative real integer scalar')
end
if ~isscalar(m) || ~all(isreal(m(:))) || ~all(mod(m(:),1)==0) || any(abs(m(:))>n)
    error('m should be a scalar integer in the range -n:n')
end
if nargin < 5 % 4th argument is really beta
    cosbeta = cos(sinbeta);
    sinbeta = sin(sinbeta);
end

% Compute legendre function for all m in 0<=m<=n:
Pn = legendre(n,cosbeta(:));

% Legendre function for |m|:
Pnm = reshape(Pn(abs(m)+1,:), size(cosbeta));

% Compute scaling term, including sign factor:
ScalingTerm = sqrt((2*n+1)/(4*pi*prod(n+(1-abs(m):abs(m))))) * (-1).^m;

% Compute exponential term:
ExpTerm = exp(1i*m*alpha);

% Put it all together:
Y = ScalingTerm .* bsxfun(@times, ExpTerm, Pnm);

% Derivatives:
if nargout>1
    if n==0
        dPmn_dbeta = zeros(size(cosbeta));
    elseif m==0
        dPmn_dbeta = Pn(2,:).';
    elseif abs(m)<n
        dPmn_dbeta = 0.5*Pn(abs(m)+2,:).' - 0.5*(n+abs(m)).*(n-abs(m)+1).*Pn(abs(m),:).';
    elseif (abs(m)==1) && (n==1)
        dPmn_dbeta = -cosbeta;
    else
        dPmn_dbeta = -abs(m).*cosbeta.*Pnm./sinbeta - (n+abs(m)).*(n-abs(m)+1).*Pn(abs(m),:).';
        dPmn_dbeta(sinbeta<=eps) = 0;
    end
    dPmn_dbeta = reshape(dPmn_dbeta, size(cosbeta));
    dY_dbeta = ScalingTerm .* bsxfun(@times, ExpTerm, dPmn_dbeta);
    dY_dalpha = Y .* 1i .* m;
end
