function i = sub2indSH(m,n)
%
% i = sub2indSH(m,n)
%
% Convert Spherical Harmonic (m,n) indices to array index i
% Assumes that i iterates from 1 (Matlab style)
%
% This file is part of the supporting data produced to accompany the
% journal article "A framework for auralization of boundary element method
% simulations including source and receiver directivity� by Jonathan A.
% Hargreaves, Luke R. Rendell and Yiu W. Lam, which was submitted for
% publication in the Journal of the Acoustical Society of America on 13th
% August 2018.
%
% The work was supported by the UK Engineering and Physical Sciences
% Research Council [grant numbers EP/J022071/1 and EP/K000012/1 �Enhanced
% Acoustic Modelling for Auralisation using Hybrid Boundary Integral
% Methods�]. It is provided without warrantee under a Creative Commons
% �Attribution� licence (CC BY); see
% http://creativecommons.org/licenses/by/4.0/ for more information.

i = n.^2 + n + m + 1;