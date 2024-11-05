function [j, djdz] = SphericalBesselReg(n,z)
%
% [j, djdz] = SphericalBesselReg(n,z)
%
% Computes a spherical Bessel function of the first kind j_n (sometimes
% called the 'regular' spherical Bessel function) and its first derivative
% if required.
%
% n must be a real scalar integer
% Returned quantities have size matching z
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
if ~isscalar(n) || ~all(isreal(n(:))) || ~all(mod(n(:),1)==0)
    error('n must be a real scalar integer')
end

% Compute j (http://dlmf.nist.gov/10.47.E3):
j = zeros(size(z));
j(abs(z(:))>eps) = sqrt(pi./(2*z(abs(z(:))>eps))) .* besselj(n+0.5,z(abs(z(:))>eps));
if any(abs(z(:))<=eps)
    if n==0
        j(abs(z)<=eps) = 1;
    else
        j(abs(z)<=eps) = 0;
    end
end

% Compute dj/dz if required (uses http://dlmf.nist.gov/10.51):
if nargout > 1
    if n==0
        djdz = -SphericalBesselReg(1,z); % This will handle the abs(z)<=eps case
    else
        djdz = zeros(size(z));
        djdz(abs(z(:))>eps) = SphericalBesselReg(n-1,z(abs(z(:))>eps)) - j(abs(z(:))>eps).*(n+1)./z(abs(z(:))>eps);
        if any(abs(z(:))<=eps)
            if n==1
                djdz(abs(z)<=eps) = 1/3;
            else
                djdz(abs(z)<=eps) = 0;
            end
        end
    end
end
