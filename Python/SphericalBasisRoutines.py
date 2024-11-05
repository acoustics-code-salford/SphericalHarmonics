# Library of spherical basus function routines containing the following functions:
#
# i = sub2indSH(m,n)
#
# (m, n) = ind2subSH(i)
#
# (r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)
#
# [nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
#
# (j, djdz) = SphericalBesselReg(n, z)
#
# j = SphericalBesselReg_pOnly(n, z)
#
# (h, dhdz) = SphericalHankelOut(n, z)
#
# h = SphericalHankelOut_pOnly(n, z)
#
# (h, dhdz) = SphericalHankelIn(n, z)
#
# h = SphericalHankelIn_pOnly(n, z)
#
# (Y, dY_dbeta, dY_dalpha) = SphericalHarmonicAll(MaxOrder, alpha, sinbeta, cosbeta)
#
# Y = SphericalHarmonicAll_pOnly(MaxOrder, alpha, sinbeta, cosbeta)
#
# (Y, dY_dbeta, dY_dalpha) = SphericalHarmonic(n, m, alpha, sinbeta, cosbeta)
#
# Y = SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta)
#
# (Phi, dPhi_dn) = SphericalBasisIn(n, m, k, x, nUV)
#
# Phi = SphericalBasisIn_pOnly(n, m, k, x)
#
# (Phi, dPhi_dn) = SphericalBasisOut(n, m, k, x, nUV)
#
# Phi = SphericalBasisIn_pOnly(n, m, k, x)
#
# (Phi, dPhi_dn) = SphericalBasisReg(n, m, k, x, nUV)
#
# Phi = SphericalBasisReg_pOnly(n, m, k, x)
#
# (Phi, dPhi_dn) = SphericalBasisOutAll(k, Bnm, x, nUV)
#
# Phi = SphericalBasisOutAll_pOnly(k, Bnm, x)
#
# Phi = pIncLoudspeaker_pOnly(LS, k, Bnm, x)
#
# Bnm = ReflectSH(Bnm, xFlag, yFlag, zFlag)
#
#
# This file is part of the supporting data produced to accompany the
# journal article "A framework for auralization of boundary element method
# simulations including source and receiver directivity” by Jonathan A.
# Hargreaves, Luke R. Rendell and Yiu W. Lam, which was submitted for
# publication in the Journal of the Acoustical Society of America on 13th
# August 2018.
#
# The work was supported by the UK Engineering and Physical Sciences
# Research Council [grant numbers EP/J022071/1 and EP/K000012/1 “Enhanced
# Acoustic Modelling for Auralisation using Hybrid Boundary Integral
# Methods”]. It is provided without warrantee under a Creative Commons
# ‘Attribution’ licence (CC BY); see
# http://creativecommons.org/licenses/by/4.0/ for more information.

from scipy.special import lpmv, spherical_jn, spherical_yn
import numpy as np

def sub2indSH(m,n):
	#
	# i = sub2indSH(m,n)
	#
	# Convert Spherical Harmonic (m,n) indices to array index i
	# Assumes that i iterates from 0 (Python style)
	i = n**2 + n + m
	return i

def ind2subSH(i):
	#
	# (m,n) = ind2subSH(i)
	#
	# Convert array index i to Spherical Harmonic (m,n) indices
	# Assumes that i iterates from 0 (Python style)
	# Assumes that arguments are NumPy arrays
	n = np.ceil(np.sqrt(i+1)-1);
	m = i - n**2 - n;	
	return (m,n)

	
def Cart2Sph(x,y,z):
	#
	# (r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)
	#
	# Converts Cartesian coordinates (x,y,z) into Spherical Polar Coordinates
	# (r, alpha, beta), where alpha is azimuth angle (angle in radians from the
	# positive x axis, with rotation around the positive z axis according to
	# the right-hand screw rule) and beta is polar angle (angle in radians from
	# the positive z axis). beta can alternatively be returned as two arrays of
	# its cos and sin values.
	#
	# It is assumed that x, y and z are all the same size.
	# The returned arrays will be the same size as the arguments.
    r = np.sqrt(x**2 + y**2 + z**2)
    rho = np.sqrt(x**2 + y**2)
    alpha = np.arctan2(y, x)
    cosbeta = z / r
    sinbeta = rho / r
	return r, alpha, sinbeta, cosbeta
	
def Cart2SphUV(x,y,z,nUV):
	#
	# [nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
	#
	# Returns the dot product between a given unit vector nUV and the unit 
	# vectors of the spherical polar coordinate system (r, alpha, beta), 
	# evaluated at Cartesian coordinates (x,y,z). Here alpha is
	# azimuth angle (angle in radians from the positive x axis, with rotation
	# around the positive z axis according to the right-hand screw rule) and
	# beta is polar angle (angle in radians from the positive z axis).
	#
	# It is assumed that x, y and z are all the same size.
	# The returned arrays will be the same size as the arguments.
    r = np.sqrt(x**2 + y**2 + z**2)
    rho = np.sqrt(x**2 + y**2)
    nUVrUV = nUV[0] * x / r + nUV[1] * y / r + nUV[2] * z / r
    nUValphaUV = nUV[1] * x / rho - nUV[0] * y / rho
    nUVbetaUV = nUV[0] * xz * x / (rho * r) + nUV[1] * z * y / (rho * r) - nUV[2] * rho / r
	return nUVrUV, nUValphaUV, nUVbetaUV
	
def SphericalBesselReg(n, z):
	#
	# (j, djdz) = SphericalBesselReg(n, z)
	# 
	# Computes a spherical Bessel function of the first kind j_n (sometimes
	# called the 'regular' spherical Bessel function) and its first derivative.
    j = spherical_jn(n,z,False)
    djdz = spherical_jn(n,z,True)
    return (j, djdz)

def SphericalHankelOut(n, z):
	#
	# (h, dhdz) = SphericalHankelOut(n, z)
	# 
	# Computes a spherical Hankel function of the first kind (outgoing in this
	# paper's lingo) and its first derivative.
    h = spherical_jn(n,z,False) + 1j*spherical_yn(n,z,False)
    dhdz = spherical_jn(n,z,True) + 1j*spherical_yn(n,z,True)
    return (h, dhdz)

def SphericalHankelIn(n, z):
	#
	# (h, dhdz) = SphericalHankelIn(n, z)
	# 
	# Computes a spherical Hankel function of the second kind (incoming in this
	# paper's lingo) and its first derivative.
    h = spherical_jn(n,z,False) - 1j*spherical_yn(n,z,False)
    dhdz = spherical_jn(n,z,True) - 1j*spherical_yn(n,z,True)
    return (h, dhdz)

def SphericalBesselReg_pOnly(n, z):
	#
	# j = SphericalBesselReg_pOnly(n, z)
	# 
	# Computes a spherical Bessel function of the first kind j_n 
	# (sometimes called the 'regular' spherical Bessel function).
	# Identical to SphericalBesselReg but does not compute the derivative.

    h = spherical_jn(n,z,False)
    return h

def SphericalHankelOut_pOnly(n, z):
	#
	# h = SphericalHankelOut_pOnly(n, z)
	# 
	# Computes a spherical Hankel function of the first kind
	# (outgoing in this paper's lingo).
	# Identical to SphericalHankelOut but does not compute the derivative.
    h = spherical_jn(n,z,False) + 1j*spherical_yn(n,z,False)
    return h

def SphericalHankelIn_pOnly(n, z):
	#
	# h = SphericalHankelIn_pOnly(n, z)
	# 
	# Computes a spherical Hankel function of the second kind 
	# (incoming in this paper's lingo).
	# Identical to SphericalHankelIn but does not compute the derivative.
    h = spherical_jn(n,z,False) - 1j*spherical_yn(n,z,False)
    return h

def SphericalHarmonicAll(MaxOrder, alpha, sinbeta, cosbeta):
	#
	# (Y, dY_dbeta, dY_dalpha) = SphericalHarmonicAll(MaxOrder, alpha, sinbeta, cosbeta)
	#
	# Computes a Spherical Harmonic function and it's angular derivatives for
	# all (m,n) up to the given maximum order. The algorithm is equivalent to that
	# implemented in SphericalHarmonic, but this version avoids repeated calls
	# to lpmv, since that is very time consuming.
	#
	# Arguments - these should all be scalars:
	# r is radius
	# alpha is azimuth angle (angle in radians from the positive x axis, with
	# rotation around the positive z axis according to the right-hand screw rule)
	# beta is polar angle, but it is specified as two arrays of its cos and sin values. 
	# MaxOrder is maximum Spherical Harmonic order and should be a non-negative real integer scalar
	#
	# Returned data will be vectors of length (MaxOrder+1)^2.
    
    # Preallocate output arrays:
    Y =  np.zeros((MaxOrder+1)**2, np.complex128)
    dY_dbeta = np.zeros((MaxOrder+1)**2, np.complex128)
    dY_dalpha = np.zeros((MaxOrder+1)**2, np.complex128)
    
    #% Loop over n and calculate spherical harmonic functions Ynm
    for n in range(MaxOrder+1):
    
        # Compute Legendre function and its derivatives for all m:
        Pn = lpmv(range(0,n+1), n, cosbeta)
        
        for m in range(-n, n+1):
    
            # Legendre function its derivatives for |m|:
            Pnm = Pn[abs(m)]
            if n==0:
                dPmn_dbeta = 0
            elif m==0:
                dPmn_dbeta = Pn[1]
            elif abs(m)<n:
                dPmn_dbeta = 0.5*Pn[abs(m)+1] - 0.5*(n+abs(m))*(n-abs(m)+1)*Pn[abs(m)-1];
            elif (abs(m)==1) and (n==1):
                dPmn_dbeta = -cosbeta
            elif sinbeta<=np.finfo(float).eps:
                dPmn_dbeta = 0
            else:
                dPmn_dbeta = -abs(m)*cosbeta*Pnm/sinbeta - (n+abs(m))*(n-abs(m)+1)*Pn[abs(m)-1]
    
            # Compute scaling term, including sign factor:
            ScalingTerm = ((-1)**m) * np.sqrt((2 * n + 1) / (4 * np.pi * np.prod(np.float64(range(n-abs(m)+1, n+abs(m)+1)))))
    
            # Compute exponential term:
            ExpTerm = np.exp(1j*m*alpha)
            
            # Put it all together:
            i = sub2indSH(m,n)
            Y[i] = ScalingTerm * ExpTerm * Pnm
            dY_dbeta[i] = ScalingTerm * ExpTerm * dPmn_dbeta
            dY_dalpha[i] = Y[i] * 1j * m
            
    return (Y, dY_dbeta, dY_dalpha)


def SphericalHarmonicAll_pOnly(MaxOrder, alpha, sinbeta, cosbeta):
	#
	# Y = SphericalHarmonicAll_pOnly(MaxOrder, alpha, sinbeta, cosbeta)
	#
	# Computes a Spherical Harmonic function for all (m,n) up to the given maximum order.
	# Identical to SphericalHarmonicAll but does not compute the derivatives.
	# The algorithm is equivalent to that implemented in SphericalHarmonic, but this
	# version avoids repeated calls to lpmv, since that is very time consuming.
	#
	# Arguments - these should all be scalars:
	# r is radius
	# alpha is azimuth angle (angle in radians from the positive x axis, with
	# rotation around the positive z axis according to the right-hand screw rule)
	# beta is polar angle, but it is specified as two arrays of its cos and sin values. 
	# MaxOrder is maximum Spherical Harmonic order and should be a non-negative real integer scalar
	#
	# Returned data will be a vector of length (MaxOrder+1)^2.

    # Preallocate output arrays:
    Y = np.zeros((MaxOrder + 1) ** 2, np.complex128)

    # % Loop over n and calculate spherical harmonic functions Ynm
    for n in range(MaxOrder + 1):

        # Compute associated Legendre function and its derivatives for all m:
        Pn = lpmv(range(0, n + 1), n, cosbeta)

        for m in range(-n, n + 1):

            # Associated Legendre function for |m|:
            Pnm = Pn[abs(m)]

            # Compute scaling term, including sign factor:
            ScalingTerm = ((-1) ** m) * np.sqrt((2 * n + 1) / (4 * np.pi * np.prod(np.float64(range(n - abs(m) + 1, n + abs(m) + 1)))))

            # Compute exponential term:
            ExpTerm = np.exp(1j * m * alpha)

            # Put it all together:
            i = sub2indSH(m,n)
            Y[i] = ScalingTerm * ExpTerm * Pnm

    return Y


def SphericalHarmonic(n, m, alpha, sinbeta, cosbeta):
	#
	# (Y, dY_dbeta, dY_dalpha) = SphericalHarmonic(n, m, alpha, sinbeta, cosbeta)
	#
	# Computes a Spherical Harmonic function of order (m,n) and it's angular derivatives.
	#
	# Arguments - these should all be scalars:
	# r is radius
	# alpha is azimuth angle (angle in radians from the positive x axis, with
	# rotation around the positive z axis according to the right-hand screw rule)
	# beta is polar angle, but it is specified as two arrays of its cos and sin values. 
	# m and n should be integer scalars; n should be non-negative and m should be in the range -n<=m<=n
	#
	# Returned data will be vectors of length (Order+1)^2.

    # Associated Legendre function its derivatives for |m|:
    Pnm = lpmv(abs(m), n, cosbeta)
    if n == 0:
        dPmn_dbeta = 0
    elif m == 0:
        dPmn_dbeta = lpmv(1, n, cosbeta)
    elif abs(m) < n:
        dPmn_dbeta = 0.5 * lpmv(abs(m) + 1, n, cosbeta) - 0.5 * (n + abs(m)) * (n - abs(m) + 1) * lpmv(abs(m) - 1, n, cosbeta);
    elif (abs(m) == 1) and (n == 1):
        dPmn_dbeta = -cosbeta
    elif sinbeta<=np.finfo(float).eps:
        dPmn_dbeta = 0
    else:
        dPmn_dbeta = -abs(m) * cosbeta * Pnm / sinbeta - (n + abs(m)) * (n - abs(m) + 1) * lpmv(abs(m) - 1, n, cosbeta)

    # Compute scaling term, including sign factor:
    ScalingTerm = ((-1) ** m) * np.sqrt((2 * n + 1) / (4 * np.pi * np.prod(np.float64(range(n - abs(m) + 1, n + abs(m) + 1)))))

    # Compute exponential term:
    ExpTerm = np.exp(1j * m * alpha)

    # Put it all together:
    Y = ScalingTerm * ExpTerm * Pnm
    dY_dbeta = ScalingTerm * ExpTerm * dPmn_dbeta
    dY_dalpha = Y * 1j * m

    return (Y, dY_dbeta, dY_dalpha)


def SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta):
	#
	# Y = SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta)
	#
	# Computes a Spherical Harmonic function of order (m,n).
	# Identical to SphericalHarmonic but does not compute the derivatives.
	#
	# Arguments - these should all be scalars:
	# r is radius
	# alpha is azimuth angle (angle in radians from the positive x axis, with
	# rotation around the positive z axis according to the right-hand screw rule)
	# beta is polar angle, but it is specified as two arrays of its cos and sin values. 
	# m and n should be integer scalars; n should be non-negative and m should be in the range -n<=m<=n
	#
	# Returned data will be a vector of length (Order+1)^2.

    # Associated Legendre function its derivatives for |m|:
    Pnm = lpmv(abs(m), n, cosbeta)

    # Compute scaling term, including sign factor:
    ScalingTerm = ((-1)**m) * np.sqrt((2 * n + 1) / (4 * np.pi * np.prod(np.float64(range(n-abs(m)+1, n+abs(m)+1)))))

    # Compute exponential term:
    ExpTerm = np.exp(1j*m*alpha)
    
    # Put it all together:
    Y = ScalingTerm * ExpTerm * Pnm

    return Y


def SphericalBasisIn(n, m, k, x, nUV):
	#
	# (Phi, dPhi_dn) = SphericalBasisIn(n, m, k, x, nUV)
	#
	# Returns Phi and dPhi/dn for an incoming Spherical Basis function of order (m,n).
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	# nUV   unit vector defining direction in which to compute dPhi/dn - 1x3
	#
	# Returned quantities are vectors with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # dot products of nUV with unit vectors of spherical coordinate system (at x):
	[nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
	
    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    (Y, dY_dbeta, dY_dalpha) = SphericalHarmonic(n, m, alpha, sinbeta, cosbeta)
    (R, dR_dkr) = SphericalHankelIn(n, k*r)

	# Evaluate Phi and dPhi/dn:
    Phi = R * Y
    dPhi_dn = (nUVrUV * k * dR_dkr * Y + (R / r) * (nUVbetaUV * dY_dbeta + nUValphaUV * dY_dalpha / sinbeta))
	
    return (Phi, dPhi_dn)


def SphericalBasisIn_pOnly(n, m, k, x):
	#
	# Phi = SphericalBasisIn_pOnly(n, m, k, x)
	#
	# Returns Phi for an incoming Spherical Basis function of order (m,n).
	# Identical to SphericalBasisIn but does not compute the derivatives.
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	#
	# Returned quantity is a vector with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)
	
    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    Y = SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta)
    R = SphericalHankelIn_pOnly(n, k*r)

	# Evaluate Phi:
    Phi = R * Y
	
    return Phi


def SphericalBasisOut(n, m, k, x, nUV):
	#
	# (Phi, dPhi_dn) = SphericalBasisOut(n, m, k, x, nUV)
	#
	# Returns Phi and dPhi/dn for an outgoing Spherical Basis function of order (m,n).
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	# nUV   unit vector defining direction in which to compute dPhi/dn - 1x3
	#
	# Returned quantities are vectors with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # dot products of nUV with unit vectors of spherical coordinate system (at x):
	[nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
	
    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    (Y, dY_dbeta, dY_dalpha) = SphericalHarmonic(n, m, alpha, sinbeta, cosbeta)
    (R, dR_dkr) = SphericalHankelOut(n, k*r)

	# Evaluate Phi and dPhi/dn:
    Phi = R * Y
    dPhi_dn = (nUVrUV * k * dR_dkr * Y + (R / r) * (nUVbetaUV * dY_dbeta + nUValphaUV * dY_dalpha / sinbeta))
	
    return (Phi, dPhi_dn)


def SphericalBasisOut_pOnly(n, m, k, x):
	#
	# Phi = SphericalBasisIn_pOnly(n, m, k, x)
	#
	# Returns Phi for an outgoing Spherical Basis function of order (m,n).
	# Identical to SphericalBasisOut but does not compute the derivatives.
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	#
	# Returned quantity is a vector with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    Y = SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta)
    R = SphericalHankelOut_pOnly(n, k*r)

	# Evaluate Phi:
    Phi = R * Y
	
    return Phi


def SphericalBasisReg(n, m, k, x, nUV):
	#
	# (Phi, dPhi_dn) = SphericalBasisReg(n, m, k, x, nUV)
	#
	# Returns Phi and dPhi/dn for an regular Spherical Basis function of order (m,n).
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	# nUV   unit vector defining direction in which to compute dPhi/dn - 1x3
	#
	# Returned quantities are vectors with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # dot products of nUV with unit vectors of spherical coordinate system (at x):
	[nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
	
    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    (Y, dY_dbeta, dY_dalpha) = SphericalHarmonic(n, m, alpha, sinbeta, cosbeta)
    (R, dR_dkr) = SphericalBesselReg(n, k*r)

	# Evaluate Phi and dPhi/dn:
    Phi = R * Y
    dPhi_dn = (nUVrUV * k * dR_dkr * Y + (R / r) * (nUVbetaUV * dY_dbeta + nUValphaUV * dY_dalpha / sinbeta))
	
    return (Phi, dPhi_dn)


def SphericalBasisReg_pOnly(n, m, k, x):
	#
	# Phi = SphericalBasisReg_pOnly(n, m, k, x)
	#
	# Returns Phi for an regular Spherical Basis function of order (m,n).
	# Identical to SphericalBasisReg but does not compute the derivatives.
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# n     must be a non-negative real integer scalar
	# m     m must be a real integer scalar in the range -n to +n
	# x     evaluation positions - real-valued array with 3 columns
	#
	# Returned quantity is a vector with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # Evaluate spherical harmonic and Hankel functions and their derivatives:
    Y = SphericalHarmonic_pOnly(n, m, alpha, sinbeta, cosbeta)
    R = SphericalHankelIn_pOnly(n, k*r)

	# Evaluate Phi:
    Phi = R * Y
	
    return Phi


def SphericalBasisOutAll(k, Bnm, x, nUV):
	#
	# (Phi, dPhi_dn) = SphericalBasisOutAll(k, Bnm, x, nUV)
	#
	# Returns Phi and dPhi/dn for a summation of outgoing Spherical Basis functions with
	# coefficients Bnm. The algorithm is equivalent to that implemented in SphericalBasisOut,
	# but this version avoids repeated call to lpmv, since that is very time consuming.
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# Bnm   directivity coefficients - vector with a square number of elements
	# x     evaluation positions - real-valued array with 3 columns
	# nUV   unit vector defining direction in which to compute dPhi/dn - 1x3
	#
	# Returned quantities are vectors with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # dot products of nUV with unit vectors of spherical coordinate system (at x):
	[nUVrUV, nUValphaUV, nUVbetaUV] = Cart2SphUV(x,y,z,nUV)
	
    # Evaluate spherical harmonic functions and their derivatives:
    if (Bnm.ndim != 1):
        raise IndexError('Bnm must be 1-dimensional')
    Order = int(np.sqrt(Bnm.size)) - 1
    (Y, dY_dbeta, dY_dalpha) = SphericalHarmonicAll(Order, alpha, sinbeta, cosbeta)

    # Loop over m and n and evalute Phi and dPhi/dn:
    Phi = np.zeros(r.size, np.complex128)
    dPhi_dn = np.zeros(r.size, np.complex128)
    for n in range(Order + 1):
        (R, dR_dkr) = SphericalHankelOut(n, k*r)
        for m in range(-n, n + 1):
            i = sub2indSH(m,n)
            Phi += Bnm[i] * R * Y[i]
            dPhi_dn += Bnm[i] * (nUVrUV * k * dR_dkr * Y[i] + (R / r) * (nUVbetaUV * dY_dbeta[i] + nUValphaUV * dY_dalpha[i] / sinbeta))

    return (Phi, dPhi_dn)	


def SphericalBasisOutAll_pOnly(k, Bnm, x, nUV):
	#
	# Phi = SphericalBasisOutAll_pOnly(k, Bnm, x)
	#
	# Returns Phi for a summation of outgoing Spherical Basis functions with coefficients Bnm.
	# The algorithm is equivalent to that implemented in SphericalBasisOut_pOnly,
	# but this version avoids repeated call to lpmv, since that is very time consuming.
	# Identical to SphericalBasisOutAll but does not compute the derivatives.
	#
	# Arguments:
	# k     wavenumber - positive real scalar
	# Bnm   directivity coefficients - vector with a square number of elements
	# x     evaluation positions - real-valued array with 3 columns
	#
	# Returned quantities is a vector with the same number of elements as x has rows.

    # Convert cartesison coordinates x to spherical coordinates:
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(x,y,z)

    # Evaluate spherical harmonic functions and their derivatives:
    if (Bnm.ndim != 1):
        raise IndexError('Bnm must be 1-dimensional')
    Order = int(np.sqrt(Bnm.size)) - 1
    Y = SphericalHarmonicAll_pOnly(Order, alpha, sinbeta, cosbeta)

    # Loop over m and n and evalute Phi:
    Phi = np.zeros(r.size, np.complex128)
    for n in range(Order + 1):
        R = SphericalHankelOut_pOnly(n, k*r)
        for m in range(-n, n + 1):
            i = sub2indSH(m,n)
            Phi += Bnm[i] * R * Y[i]

    return Phi	
	
	

def pIncLoudspeaker_pOnly(LS, k, Bnm, x):
	#
	# Phi = pIncLoudspeaker_pOnly(LS, k, Bnm, x)
	#
	# Returns Phi for a summation of outgoing Spherical Basis functions with coefficients Bnm.
	# Very similar to SphericalBasisOutAll_pOnly but differs in that the loudspeaker coordinate
	# system may be rotated and is in front-pole format.
	#
	# Arguments:
	# LS 	structure with four fields 'center', 'front', 'top' and 'side', all of which are 3 element vectors.
	# k     wavenumber - positive real scalar
	# Bnm   directivity coefficients - vector with a square number of elements
	# x     evaluation positions - real-valued array with 3 columns
	#
	# Returned quantities is a vector with the same number of elements as x has rows.

    # vector from LS center to x in cartesian coordinate system rotated with loudspeaker. Note that Bnm is in front-pole format. 
	(r, alpha, sinbeta, cosbeta) = Cart2Sph(np.dot(x - LS.center, LS.top), np.dot(x - LS.center, LS.side), np.dot(x - LS.center, LS.front))

    # Evaluate spherical harmonic functions and their derivatives:
    if (Bnm.ndim != 1):
        raise IndexError('Bnm must be 1-dimensional')
    Order = int(np.sqrt(Bnm.size)) - 1
    Y = SphericalHarmonicAll_pOnly(Order, alpha, sinbeta, cosbeta)

    # Loop over m and n and evalute Phi:
    Phi = np.zeros(r.size, np.complex128)
    for n in range(Order + 1):
        R = SphericalHankelOut_pOnly(n, k*r)
        for m in range(-n, n + 1):
            i = sub2indSH(m,n)
            Phi += Bnm[i] * R * Y[i]

    return Phi



def ReflectSH(Bnm, xFlag, yFlag, zFlag):
	#
	# Bnm = ReflectSH(Bnm, xFlag, yFlag, zFlag)   
	#
	# Reflect an Spherical Harmonic representation of a sound-field in 1 to 3 cartesian axes.
	#
	# Argumments:
	# Bnm       Vector of Spherical Harmonic weights. Must have (Order+1)^2 entries, where Order is an integer.
	# xFlag		Boolean indicating whether to flip in the x-direction
	# yFlag		Boolean indicating whether to flip in the y-direction
	# zFlag		Boolean indicating whether to flip in the z-direction

	# Get lists of n and m values:
	(m, n) = ind2subSH(np.arange(Bnm.size))

	# Reflecting in Z:
	if zFlag:
		Bnm = Bnm * (-1.)**(n+m)

	# Reflecting in X:
	if xFlag:
		Bnm = Bnm * (-1.)**m

	# Reflecting in X or Y:
	if xFlag ^ yFlag # XOR
		for n in range(int(ceil(sqrt(Bnm.size)))-1)
		   i = sub2indSH(np.arange(-n,n+1),n)
		   Bnm[i] = np.flip(Bnm[i])
		end
	end
	
	return Bnm
	