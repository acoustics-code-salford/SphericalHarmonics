function T = GetTranslationMatrix(t,k,OrderS,OrderR,WaveTypes)
%
% T = GetTranslationMatrix(t,k,OrderS,OrderR,WaveType)
%
% Computes a translation matrix T from the coefficients of a Spherical
% Harmonic source (outgoing spherical Hankel radial functions) to the
% coefficients at a Spherical Harmonic receiver (spherical Bessel radial
% functions) location at position t relative to the source. It is assumed
% that both spherical coordinate systems (source and receiver) are aligned
% to the same Cartesian system in which t is expressed. a is the polar
% angle from the postive z axis.
%
% Essentially computes equation 3.2.17 of: 
% Gumerov, N., & Duraiswami, R. (2005). Fast Multipole Methods for the
% Helmholtz Equation in Three Dimensions (1st ed.). Elsevier Science.
%
% Modified 3rd Nov 2023 to allow differing wave types, essentially adding
% eq. 3.2.18 to the mix
%
% Arguments:
%
% t         Cartesian translation vector (1x3 real row vector)
% k         Wavenumber (positive real scalar or vector in radians/meter)
% OrderS    Order of the source (non-negative real integer scalar)
% OrderR    Order of the receiver (non-negative real integer scalar)
% WaveTypes Should  be equal to 'Out2Out', 'Out2Reg', or 'Reg2Reg'
%
% This file also contains the sub-functions
% GetStructuralTranslationCoefficients and Wigner3jSymbol.
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
    if ~isrow(t) || (numel(t)~=3) || any(imag(t(:))~=0)
        error('t must be a real 1x3 row vector')
    end
    if ~isvector(k) || any(k(:)<=0) || any(imag(k(:))~=0)
        error('k must be a positive real scalar or vector')
    end
    if ~isscalar(OrderS) || any(OrderS(:)<0) || any(imag(OrderS(:))~=0) || any(mod(OrderS(:),1)~=0)
        error('OrderS must be a non-negative real integer scalar')
    end
    if ~isscalar(OrderR) || any(OrderR(:)<0) || any(imag(OrderR(:))~=0) || any(mod(OrderR(:),1)~=0)
        error('OrderR must be a non-negative real integer scalar')
    end
    OrderT = OrderS + OrderR;
    if nargin < 5
        WaveTypes = 'Out2Reg'; % the default in the published version
    elseif ~any(strcmp(WaveTypes,{'Out2Out','Out2Reg','Reg2Reg'}))
        error('WaveTypes must be equal to ''Out2Out'', ''Out2Reg'', or ''Reg2Reg''');
    end

    % Pre-calculate table of structural coefficients:
    % Variable is persistent and stored between calls to avoid recalculation.
    % It is only recalculated if the values of OrderS and OrderR change.
    persistent S
    if isempty(S)
        S = GetStructuralTranslationCoefficients(OrderS,OrderR);
    elseif (size(S,1)~=(OrderR+1)^2) || (size(S,2)~=(OrderS+1)^2) || (size(S,3)~=(OrderR+OrderS+1)^2)
        S = GetStructuralTranslationCoefficients(OrderS,OrderR);
    end    

    % Express t in spherical coordinates:
    [r,alpha,sinbeta,cosbeta] = Cart2Sph(t(1),t(2),t(3));

    % Evaluate spherical harmonic functions:
    Y = SphericalHarmonicAll(OrderT,alpha,sinbeta,cosbeta);

    % Allocated results array:
    T = zeros((OrderR+1)^2, (OrderS+1)^2, numel(k));

    % Loop over entries in k and translation order:
    for ik = 1:numel(k)

        % Loop over translation order & compute summation:
        for nT = 0:OrderT
            if strcmp(WaveTypes, 'Out2Reg')
                fr = SphericalHankelOut(nT,k(ik)*r); % Compute radial function:
            else
                fr = SphericalBesselReg(nT,k(ik)*r); % Compute radial function:
            end
            for mT = -nT:nT
                iT = sub2indSH(mT,nT);
                T(:,:,ik) = T(:,:,ik) + fr * Y(iT) * S(:,:,iT);
            end
        end
    end

end

function S = GetStructuralTranslationCoefficients(OrderS,OrderR)
%
% S = GetStructuralTranslationCoefficients(OrderS,OrderR)
%
% Computes the 'Structural Translation Coefficients' used in Spherical
% Harmonic translation routines, as defined in section 3.2.1 of: 
% Gumerov, N., & Duraiswami, R. (2005). Fast Multipole Methods for the
% Helmholtz Equation in Three Dimensions (1st ed.). Elsevier Science.
%
% Arguments:
% OrderS    Order of the source   (non-negative real integer scalar)
% OrderR    Order of the receiver (non-negative real integer scalar)
% 
% Returned variable is a 3D array of size [(OrderR+1)^2, (OrderS+1)^2,
% (OrderR+OrderS+1)^2]. 

    % Argument checking:
    if ~isscalar(OrderS) || any(OrderS(:)<0) || any(imag(OrderS(:))~=0) || any(mod(OrderS(:),1)~=0)
        error('OrderS must be a non-negative real integer scalar')
    end
    if ~isscalar(OrderR) || any(OrderR(:)<0) || any(imag(OrderR(:))~=0) || any(mod(OrderR(:),1)~=0)
        error('OrderR must be a non-negative real integer scalar')
    end

    % Order required for translation:
    OrderT = OrderS + OrderR;

    % Allocate cell array:
    S = zeros((OrderR+1)^2, (OrderS+1)^2, (OrderT+1)^2,1);

    % Loop over translation order (n2 & m2):
    h = waitbar(0, 'Computing Structural Translation Coefficients');
    for nT = 0:OrderT % n'' in book
        for mT = -nT:nT % m'' in book
            iT = sub2indSH(mT,nT);
            if mT < 0 % because m'' is negated
                epT = (-1)^mT;
            else
                epT = 1;
            end

            % Loop over source order (nS & mS):
            for nS = 0:OrderS % n in book
                for mS = -nS:nS % m in book
                    if mS > 0
                        epS = (-1)^mS;
                    else
                        epS = 1;
                    end

                    % Loop over recevier order (nR & mR):
                    for nR = 0:OrderR % n' in book
                        for mR = -nR:nR % m' in book
                            if mR < 0 % because m' is negated
                                epR = (-1)^mR;
                            else
                                epR = 1;
                            end

                            % Compute coefficient if within non-zero range:
                            if (nT>=abs(nR-nS)) && (nT<=(nR+nS))
                                S(sub2indSH(mR,nR), sub2indSH(mS,nS), iT) = ...
                                            1i^(nR+nT-nS) * epS * epR * epT ...
                                            * sqrt(4*pi*(2*nS+1)*(2*nR+1)*(2*nT+1))...
                                            * Wigner3jSymbol(nS, nR, nT, mS, -mR, -mT)...
                                            * Wigner3jSymbol(nS, nR, nT, 0, 0, 0);
                            end
                        end
                    end

                end
            end
            waitbar(iT/((OrderT+1)^2), h);
        end
    end
    delete(h);

end        
        
        
function W3jS = Wigner3jSymbol(j1, j2, j3, m1, m2, m3)

    % W3jS = Wigner3j(j1, j2, j3, m1, m2, m3)
    %
    % Computes the Wigner 3j symbol following the formulaetion given at
    % http://mathworld.wolfram.com/Wigner3j-Symbol.html.
    %
    % Arguments:
    %
    % j1, j2, j3, m1, m2 and m3     All must be scalar half-integers

    % Check arguments against 'selection rules' (cited to Messiah 1962, pp. 1054-1056; Shore and Menzel 1968, p. 272)
    % Nullifying any of these means the symbol equals zero.
    if (abs(m1)<=abs(j1)) && (abs(m2)<=abs(j2)) && (abs(m3)<=abs(j3)) && (m1+m2+m3==0) && (abs(j1-j2)<=j3) && (j3<=(j1+j2)) && (mod(j1+j2+j3,1)==0)

        % Evaluate the symbol using the Racah formula (Equation 7):

        % Evalaute summation:
        W3jS = 0;
        for t = 0:min([j1+j2-j3, j1-m1, j2+m2])
            if (j3-j2+t+m1>=0) && (j3-j1+t-m2>=0) && (j1+j2-j3-t>=0) && (j1-t-m1>=0) && (j2-t+m2>=0)
                % Only include term in summation if all factorials have non-negative arguments
                x = factorial(t)...
                  * factorial(j3-j2+t+m1)...
                  * factorial(j3-j1+t-m2)...
                  * factorial(j1+j2-j3-t)...
                  * factorial(j1-t-m1)...
                  * factorial(j2-t+m2);
                W3jS = W3jS + (-1)^t/x;
            end
        end

        % Coefficients outside the summation:
        W3jS = W3jS...
             * (-1)^(j1-j2-m3)...
             * sqrt(factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2)* factorial(j3+m3)*factorial(j3-m3))...
             * sqrt(factorial(j1 + j2 - j3)*factorial(j1 - j2 + j3)*factorial(-j1 + j2 + j3) / factorial(j1 + j2 + j3 + 1));

    else
        W3jS = 0; % One of the 'Selection Rules' was nullified.
    end

end
