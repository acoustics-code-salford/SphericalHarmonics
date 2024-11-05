function [rUV, alphaUV, betaUV] = Cart2SphUV(x,y,z,nUV)

% [rUV, alphaUV, betaUV] = Cart2SphUV(x,y,z)
% [rUV.nUV, alphaUV.nUV, betaUV.nUV] = Cart2SphUV(x,y,z,nUV)
%
% Returns the unit vectors of the spherical polar coordinate system (r,
% alpha, beta), evaluated at Cartesian coordinates (x,y,z). Here alpha is
% azimuth angle (angle in radians from the positive x axis, with rotation
% around the positive z axis according to the right-hand screw rule) and
% beta is polar angle (angle in radians from the positive z axis).
% Alternatively, if a unit vector nUV is given as a fourth argument, then
% the dot product between this and the unit vectors of the spherical polar
% coordinate system are returned instead.
%
% It is assumed that x, y and z are all the same size.
% The return arrays will all be size numel(x) x 3 if nUV is not
% provided, or size(x) if it is.
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

% Spherical and cylindrical radii:
r = sqrt(x(:).^2 + y(:).^2 + z(:).^2);
rho = sqrt(x(:).^2 + y(:).^2);

% Compute unit vectors:
rUV = [x(:)./r, y(:)./r, z(:)./r];
alphaUV = [-y(:)./rho, x(:)./rho, zeros(size(z(:)))];
betaUV = [x(:).*z(:)./(rho.*r), y(:).*z(:)./(rho.*r), -rho./r];
if any(rho<=eps) % Catch cases on polar axis
    alphaUV(rho<=eps,:) = repmat([-1, 0, 0], nnz(rho<=eps), 1);
    betaUV(rho<=eps,:) = [0, 1, 0] .* sign(z(rho<=eps));
end



if nargin==4
    % Compute dot product between nUV and unit vectors:
    rUV     = reshape(dot3D(rUV,     nUV), size(x));
    alphaUV = reshape(dot3D(alphaUV, nUV), size(x));
    betaUV  = reshape(dot3D(betaUV,  nUV), size(x));
end
