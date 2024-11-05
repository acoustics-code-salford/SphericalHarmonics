function Bnm = ReflectSH(Bnm,AxisFlag)
%
% Bnm = ReflectSH(Bnm,AxisFlag)   
%
% Reflect an Spherical Harmonic representation of a sound-field in 1 to 3
% cartesian planes.
%
% Argumments:
% Bnm       Column vector of Spherical Harmonic weights. Must have (Order+1)^2 entries, where Order is an integer.
% AxisFlag  3 element boolean vector indicating in which directions to flip.
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

% Parameter checking:
if ~iscolumn(Bnm) || (mod(sqrt(numel(Bnm)), 1)~=0)
    error('B must be a column vector with (Order+1)^2 entries, where Order is an integer.')
end
if any(size(AxisFlag)~=[1 3])
    error('AxisFlag must be a row vector with 3 elements')
elseif any((AxisFlag~=0) & (AxisFlag~=1))
    error('AxisFlag must contain only boolean values or zeros and ones')
end

% Get n and m:
[m, n] = ind2subSH((1:numel(Bnm))');

% Reflecting in Z:
if AxisFlag(3)
    Bnm = Bnm .* (-1).^(n+m);
end

% Reflecting in X:
if AxisFlag(1)
    Bnm = Bnm .* (-1).^m;
end

% Reflecting in X or Y:
if xor(AxisFlag(1),AxisFlag(2))
    for n = 0:max(n)
       i = sub2indSH(-n:n,n);
       Bnm(i) = flipud(Bnm(i));
    end
end

    
    
