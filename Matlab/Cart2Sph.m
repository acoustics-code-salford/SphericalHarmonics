function [r,alpha,sinbeta,cosbeta] = Cart2Sph(x,y,z)

% [r,alpha,sinbeta,cosceta] = Cart2Sph(x,y,z)
% [r,alpha,beta] = Cart2Sph(x,y,z)
%
% Converts Cartesian coordinates (x,y,z) into Spherical Polar Coordinates
% (r, alpha, beta), where alpha is azimuth angle (angle in radians from the
% positive x axis, with rotation around the positive z axis according to
% the right-hand screw rule) and beta is polar angle (angle in radians from
% the positive z axis). beta can alternatively be returned as two arrays of
% its cos and sin values.
%
% It is assumed that x, y and z are all the same size.
% The returned arrays will be the same size as the arguments.
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

r = sqrt(x.^2 + y.^2 + z.^2);
rho = sqrt(x.^2 + y.^2);
alpha = atan2(y, x);
cosbeta = z ./ r;
sinbeta = rho ./ r;

if nargout<4 % 3rd output variable should be assigned the value beta
    sinbeta = atan2(sinbeta, cosbeta);
end
