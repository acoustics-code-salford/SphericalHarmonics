function [Y, dY_dbeta, dY_dalpha] = SphericalHarmonicAll(MaxOrder,alpha,sinbeta,cosbeta)
%
% [Y, dY_dbeta, dY_dalpha] = SphericalHarmonicAll(MaxOrder,alpha,beta)
% [Y, dY_dbeta, dY_dalpha] = SphericalHarmonicAll(MaxOrder,alpha,sinbeta,cosbeta)
%
% Computes a Spherical Harmonic function, and it's angular derivatives if
% required, for all (m,n) up to the given order. The algorithm is
% equivalent to that implemented in SphericalHarmonic.m, but this version
% avoids repeated call to legendre.m, since that is very time consuming.
%
% Arguments:
% r is radius
% alpha is azimuth angle (angle in radians from the positive x axis, with
% rotation around the positive z axis according to the right-hand screw rule)
% beta is polar angle (angle in radians from the positive z axis). It can
% alternatively be specified as two arrays of its cos and sin values. 
% These should all be arrays of identical size.
% MaxOrder is maximum Spherical Harmonic order and should be a non-negative real integer scalar
%
% Returned arrays will have numel(alpha) rows and (MaxOrder+1)^2 columns.
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
if ~isscalar(MaxOrder) || ~all(isreal(MaxOrder(:))) || ~all(mod(MaxOrder(:),1)==0) || any(MaxOrder(:)<0)
    error('Order must be a non-negative real integer scalar')
end
if nargin < 4 % 3rd argument is really beta
    cosbeta = cos(sinbeta);
    sinbeta = sin(sinbeta);
end

% Ensure angle data is column vectors:
alpha = alpha(:);
sinbeta = sinbeta(:);
cosbeta = cosbeta(:);

% Preallocate output arrays:
Y         = zeros(numel(alpha), (MaxOrder+1)^2);
dY_dbeta  = zeros(numel(alpha), (MaxOrder+1)^2);
dY_dalpha = zeros(numel(alpha), (MaxOrder+1)^2);


%% Loop over n and calculate spherical harmonic functions Ynm
for n = 0:MaxOrder

    % Compute legendre function for all m:
    Pn = legendre(n,cosbeta).'; % all m
    
    for m = -n:n
            
        % Legendre function for |m|:
        Pnm = Pn(:,abs(m)+1);
      
        % Compute scaling term, including sign factor:
        ScalingTerm = sqrt((2*n+1)/(4*pi*prod(n+(1-abs(m):abs(m))))) * (-1).^m;

        % Compute exponential term:
        ExpTerm = exp(1i*m*alpha);

        % Put it all together:
        i = sub2indSH(m,n);
        Y(:,i) = ScalingTerm .* ExpTerm .* Pnm;
        
        % Derivatives:
        if nargout>1
            if n==0
                dPmn_da = zeros(size(cosbeta));
            elseif m==0
                dPmn_da = Pn(:,2);
            elseif abs(m)<n
                dPmn_da = 0.5*Pn(:,abs(m)+2) - 0.5*(n+abs(m)).*(n-abs(m)+1).*Pn(:,abs(m));
            elseif (abs(m)==1) && (n==1)
                dPmn_da = -cosbeta;
            else
                dPmn_da = -abs(m).*cosbeta.*Pnm./sinbeta - (n+abs(m)).*(n-abs(m)+1).*Pn(:,abs(m));
                dPmn_da(sinbeta<=eps) = 0;
            end
            dY_dbeta(:,i) = ScalingTerm .* ExpTerm .* dPmn_da;
            dY_dalpha(:,i) = Y(:,i) .* 1i .* m;
        end          

    end
end
