function [h, dhdz] = SphericalHankelOut(n,z)
%
% [h, dhdz] = SphericalHankelOut(n,z)
%
% Computes a spherical Hankel function of the first kind (outgoing in this
% paper's lingo) and its first derivative if required.
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

% Compute h (http://dlmf.nist.gov/10.47.E3):
h = sqrt(pi./(2*z)) .* besselh(n+0.5,1,z); % first kind

% Compute dh/dz if required (uses http://dlmf.nist.gov/10.51):
if nargout > 1
    if n==0
        dhdz = -SphericalHankelOut(1,z);
    else
        dhdz = SphericalHankelOut(n-1,z) - h.*(n+1)./z;
    end
end
