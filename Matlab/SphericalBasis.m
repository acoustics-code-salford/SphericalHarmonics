function [Phi, GradPhi] = SphericalBasis(k, n, m, x, type, nUV)
%
% Phi = SphericalBasis(k, n, m, x, type)
% [Phi, GradPhi] = SphericalBasis(k, n, m, x, type)
% [Phi, dPhi_dn] = SphericalBasis(k, n, m, x, type, nUV)
%
% Returns Phi, and Grad Phi or dPhi/dn if required for Spherical Basis
% function of order (m,n).
%
% Arguments:
% k     wavenumber - positive real scalar
% n     must be a non-negative real integer scalar
% m     m must be a real integer scalar in the range -n to +n
% x     evaluation positions - real-valued array with 3 columns
% type  'In', 'Out' or 'Reg'
% nUV   unit vector defining direction in which to compute dPhi/dn - 1x3
%
% Returned arrays have the same number of rows as x, with 1 column for Phi
% and dPhi_dn or 3 columns for GradPhi.
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
if ~isscalar(k) || any(imag(k(:))~=0) || any(k(:)<=0)
    error('k must be a positive real scalar')
end
if ~isscalar(n) || any(mod(n(:),1)~=0) || any(n(:)<0)
    error('n must be a non-negative real integer scalar')
end
if ~isscalar(m) || any(mod(m(:),1)~=0) || any(abs(m(:))>n)
    error('m must be a real integer scalar in the range -n to +n')
end
if (size(x,2)~=3) || any(imag(x(:))~=0)
    error('x must be real and have 3 columns')
end
if ~any(strcmp(type,{'In','Out','Reg'}))
    error('type must be equal to ''In'', ''Out'', or ''Reg''');
end
if nargin>5
    if ~isvector(nUV) || (numel(nUV)~=3)
        error('nUV must be a 3-element vector')
    end
end

% Convert cartesian to spherical coordinates:
[r,alpha,sinbeta,cosbeta] = Cart2Sph(x(:,1), x(:,2), x(:,3));
if (nargout>1) && (nargin>5)
    [rUV, alphaUV, betaUV] = Cart2SphUV(x(:,1), x(:,2), x(:,3), nUV);
elseif nargout>1
    [rUV, alphaUV, betaUV] = Cart2SphUV(x(:,1), x(:,2), x(:,3));
end

% Evaluate spherical harmonic functions:
if nargout>1
    [Y, dY_dBeta, dY_dalpha] = SphericalHarmonic(m,n,alpha,sinbeta,cosbeta);
else
    Y = SphericalHarmonic(m,n,alpha,sinbeta,cosbeta);
end    

% Evaluate radial function:
switch type
    case 'In'
        if nargout>1
            [R, dR_dkr] = SphericalHankelIn(n,k*r);
        else
            R = SphericalHankelIn(n,k*r);
        end
    case 'Out'
        if nargout>1
            [R, dR_dkr] = SphericalHankelOut(n,k*r);
        else
            R = SphericalHankelOut(n,k*r);
        end
    case 'Reg'
        if nargout>1
            [R, dR_dkr] = SphericalBesselReg(n,k*r);
        else
            R = SphericalBesselReg(n,k*r);
        end
end

% Evaluate Phi and GradPhi:
Phi = R.*Y;
if (nargout>1) && (nargin>5)
    GradPhi = rUV.*k.*dR_dkr.*Y + (R./r) .* (betaUV.*dY_dBeta + alphaUV.*dY_dalpha./sinbeta); % Is really dPhi/dn
    if any(sinbeta<eps) % catch cases on polar axis, where NaNs can arise from derivative w.r.t. alpha:
        i0 = sinbeta<eps; % version with last term omitted:
        GradPhi(i0) = rUV(i0).*k.*dR_dkr(i0).*Y(i0) + (R(i0)./r(i0)) .* betaUV(i0).*dY_dBeta(i0); % Is really dPhi/dn
    end
elseif nargout>1
    GradPhi = rUV.*repmat(k.*dR_dkr.*Y,1,3) + repmat(R./r,1,3) .* (betaUV.*repmat(dY_dBeta,1,3) + alphaUV.*repmat(dY_dalpha./sinbeta,1,3));
    if any(sinbeta<eps) % catch cases on polar axis, where NaNs can arise from derivative w.r.t. alpha:
        i0 = sinbeta<eps; % version with last term omitted:
        GradPhi(i0,:) = rUV(i0,:).*repmat(k.*dR_dkr(i0).*Y(i0),1,3) + repmat(R(i0)./r(i0),1,3) .* betaUV(i0,:).*repmat(dY_dBeta(i0),1,3);
    end
end
