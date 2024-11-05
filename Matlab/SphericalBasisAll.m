function [Phi, GradPhi] = SphericalBasisAll(k, Bnm, x, type, nUV)
%
% Phi = SphericalBasisAll(k, Bnm, x, type)
% [Phi, GradPhi] = SphericalBasisAll(k, Bnm, x, type)
% [Phi, dPhi_dn] = SphericalBasisAll(k, Bnm, x, type, nUV)
%
% Returns Phi, and Grad Phi or dPhi/dn if required, for a summation of
% Spherical Basis functions with coefficients Bnm. The algorithm is
% equivalent to that implemented in SphericalBasis.m, but this version
% avoids repeated call to legendre.m, since that is very time consuming.
%
% Arguments:
% k     wavenumber - positive real scalar
% Bnm   directivity coefficients - vector with a square number of elements
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
if ~isvector(Bnm) || (mod(sqrt(numel(Bnm)),1)~=0)
    error('Bnm must be a vector with a square number of elements')
end
if (size(x,2)~=3) || any(imag(x(:))~=0)
    error('x must be real and have 3 columns')
end
if ~any(strcmp(type,{'In','Out','Reg'}))
    error('type must be equal to ''In'', ''Out'', or ''Reg''');
end
if nargin>4
    if ~isvector(nUV) || (numel(nUV)~=3)
        error('nUV must be a 3-element vector')
    end
end

% Convert cartesian to spherical coordinates:
[r,alpha,sinbeta,cosbeta] = Cart2Sph(x(:,1), x(:,2), x(:,3));
if (nargout>1) && (nargin>4)
    [rUV, alphaUV, betaUV] = Cart2SphUV(x(:,1), x(:,2), x(:,3), nUV);
elseif nargout>1
    [rUV, alphaUV, betaUV] = Cart2SphUV(x(:,1), x(:,2), x(:,3));
end

% Get spherical harmonic order:
Order = sqrt(numel(Bnm))-1;

% Evaluate spherical harmonic functions:
if nargout>1
    [Y, dY_dBeta, dY_dalpha] = SphericalHarmonicAll(Order,alpha,sinbeta,cosbeta);
else
    Y = SphericalHarmonicAll(Order,alpha,sinbeta,cosbeta);
end    

% Loop over m and n and evalute output data:
Phi = zeros(numel(r),1);
if (nargout>1) && (nargin>4)
    GradPhi = zeros(numel(r),1); % Is really dPhi/dn
elseif nargout>1
    GradPhi = zeros(numel(r),3);
end
for n = 0:Order
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
                [R, dR_dkr] = SphericalHankelReg(n,k*r);
            else
                R = SphericalHankelReg(n,k*r);
            end
    end
    for m = -n:n
        i = sub2indSH(m,n);
        Phi = Phi + Bnm(i).*R.*Y(:,i);
        if (nargout>1) && (nargin>4)
            if any(sinbeta<eps) % catch cases on polar axis, where NaNs can arise from derivative w.r.t. alpha:
                i0 = sinbeta>=eps; % full statement:
                GradPhi(i0) = GradPhi(i0) + Bnm(i).*(rUV(i0).*k.*dR_dkr(i0).*Y(i0,i) + (R(i0)./r(i0)) .* (betaUV(i0).*dY_dBeta(i0,i) + alphaUV(i0).*dY_dalpha(i0,i)./sinbeta(i0))); % Is really dPhi/dn
                i0 = sinbeta<eps; % version with last term omitted:
                GradPhi(i0) = GradPhi(i0) + Bnm(i).*(rUV(i0).*k.*dR_dkr(i0).*Y(i0,i) + (R(i0)./r(i0)) .* betaUV(i0).*dY_dBeta(i0,i)); % Is really dPhi/dn
            else % full statement without indexing:
                GradPhi = GradPhi + Bnm(i).*(rUV.*k.*dR_dkr.*Y(:,i) + (R./r) .* (betaUV.*dY_dBeta(:,i) + alphaUV.*dY_dalpha(:,i)./sinbeta)); % Is really dPhi/dn
            end
        elseif nargout>1
            if any(sinbeta<eps) % catch cases on polar axis, where NaNs can arise from derivative w.r.t. alpha:
                i0 = sinbeta>=eps; % full statement:
                GradPhi(i0,:) = GradPhi(i0,:) + Bnm(i).*(rUV(i0,:).*repmat(k.*dR_dkr(i0).*Y(i0,i),1,3) + repmat(R(i0)./r(i0),1,3) .* (betaUV(i0,:).*repmat(dY_dBeta(i0,i),1,3) + alphaUV(i0,:).*repmat(dY_dalpha(i0,i)./sinbeta(i0),1,3)));
                i0 = sinbeta<eps; % version with last term omitted:
                GradPhi(i0,:) = GradPhi(i0,:) + Bnm(i).*(rUV(i0,:).*repmat(k.*dR_dkr(i0).*Y(i0,i),1,3) + repmat(R(i0)./r(i0),1,3) .* betaUV(i0,:).*repmat(dY_dBeta(i0,i),1,3));
            else % full statement without indexing:
                GradPhi = GradPhi + Bnm(i).*(rUV.*repmat(k.*dR_dkr.*Y(:,i),1,3) + repmat(R./r,1,3) .* (betaUV.*repmat(dY_dBeta(:,i),1,3) + alphaUV.*repmat(dY_dalpha(:,i)./sinbeta,1,3)));
            end
        end
    end
end
