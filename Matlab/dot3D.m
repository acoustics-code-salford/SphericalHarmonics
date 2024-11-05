function c = dot3D(a,b)
%
% c = dot3D(a,b)
% 
% Returns the row-wise dot product of the matrices A and B.  
% A and B must have 3 columns. In addition, they must have an equal number
% of rows, or one of them must have just one row.
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

if (size(a,2)~=3) || (size(b,2)~=3)
    error('a and b must both have 3 columns')
elseif (size(a,1)==size(b,1))
    c = sum(a.*b, 2);
elseif isrow(a) || isrow(b)
    c = sum(bsxfun(@times,a,b), 2);
else
    error('a and b must have the same number of rows, unless one of them has only one row')
end
