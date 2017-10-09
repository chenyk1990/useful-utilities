#!/usr/bin/env python
"""
The Python Seismic Anisotropy Toolkit

:Copyright:
    Author: Lili Feng
    Graduate Research Assistant
    CIEI, Department of Physics, University of Colorado Boulder
    email: lili.feng@colorado.edu
    
:Dependencies:
    pyasdf 0.1.4
    mayavi 4.5.0
    cartopy 0.15.1
    
:References:
    Bond, W.L., 1943. The mathematics of the physical properties of crystals.
        Bell Labs Technical Journal, 22(1), pp.1-72.
    Babuska, V. and Cara, M., 1991. Seismic anisotropy in the Earth (Vol. 10).
        Springer Science & Business Media.
    Carcione, J.M., 2014. Wave fields in real media:
        Wave propagation in anisotropic, anelastic, porous and electromagnetic media (Vol. 38). Elsevier.
    Riley, K.F., Hobson, M.P. and Bence, S.J., 2006. Mathematical methods for physics and engineering: a comprehensive guide.
        Cambridge university press.
    Tsvankin, I., 2012. Seismic signatures and analysis of reflection data in anisotropic media.
        Society of Exploration Geophysicists.
    Wolfe, J.P., 2005. Imaging phonons: acoustic wave propagation in solids.
        Cambridge University Press.
::: Note :::
For rotation matrix ambiguities:
    https://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities

Direction of rotation:
    The direction of rotation is given by the right-hand rule (orient the thumb of the right hand along the axis around which the rotation occurs,
    with the end of the thumb at the positive end of the axis; curl your fingers; the direction your fingers curl is the direction of rotation).
    Therefore, the rotations are counterclockwise if looking along the axis of rotation from positive to negative.
"""
import numpy as np
import copy
import pyasdf
try:
    from opt_einsum import contract
    use_opt_einsum=True
except: use_opt_einsum=False
import mayavi.mlab
import cartopy.crs as ccrs
from matplotlib import cm
import matplotlib.pyplot as plt
from numba import jit, float32, int32

########################################################################################################
# axis sequences for Euler angles
_NEXT_AXIS = [1, 2, 0, 1]

# map axes strings to/from tuples of inner axis, parity, repetition, frame
_AXES2TUPLE = {
    'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
    'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
    'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
    'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
    'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
    'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
    'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
    'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

_TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())

# For testing whether a number is close to zero
_EPS4 = np.finfo(float).eps * 4.0
########################################################################################################


def euler2mat(ai, aj, ak, axes='sxyz'):
    """Return rotation matrix from Euler angles and axis sequence.
    Parameters
    ----------
    ai : float
        First rotation angle in degree (according to `axes`).
    aj : float
        Second rotation angle in degree (according to `axes`).
    ak : float
        Third rotation angle in degree (according to `axes`).
    axes : str, optional
        Axis specification; one of 24 axis sequences as string or encoded
        tuple - e.g. ``sxyz`` (the default).
    Returns
    -------
    mat : array-like shape (3, 3) or (4, 4)
        Rotation matrix or affine.
    Examples
    --------
    >>> R = euler2mat(1, 2, 3, 'syxz')
    >>> np.allclose(np.sum(R[0]), -1.34786452)
    True
    >>> R = euler2mat(1, 2, 3, (0, 1, 0, 1))
    >>> np.allclose(np.sum(R[0]), -0.383436184)
    True
    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes
    ai=np.pi*ai/180.; aj=np.pi*aj/180.; ak=np.pi*ak/180.
    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = np.sin(ai), np.sin(aj), np.sin(ak)
    ci, cj, ck = np.cos(ai), np.cos(aj), np.cos(ak)
    cc, cs = ci*ck, ci*sk
    sc, ss = si*ck, si*sk

    M = np.eye(3)
    if repetition:
        M[i, i] = cj
        M[i, j] = sj*si
        M[i, k] = sj*ci
        M[j, i] = sj*sk
        M[j, j] = -cj*ss+cc
        M[j, k] = -cj*cs-sc
        M[k, i] = -sj*ck
        M[k, j] = cj*sc+cs
        M[k, k] = cj*cc-ss
    else:
        M[i, i] = cj*ck
        M[i, j] = sj*sc-cs
        M[i, k] = sj*cc+ss
        M[j, i] = cj*sk
        M[j, j] = sj*ss+cc
        M[j, k] = sj*cs-sc
        M[k, i] = -sj
        M[k, j] = cj*si
        M[k, k] = cj*ci
    return M

def mat2euler(mat, axes='sxyz'):
    """Return Euler angles from rotation matrix for specified axis sequence.
    Note that many Euler angle triplets can describe one matrix.
    Parameters
    ----------
    mat : array-like shape (3, 3) or (4, 4)
        Rotation matrix or affine.
    axes : str, optional
        Axis specification; one of 24 axis sequences as string or encoded
        tuple - e.g. ``sxyz`` (the default).
    Returns
    -------
    ai : float
        First rotation angle (according to `axes`).
    aj : float
        Second rotation angle (according to `axes`).
    ak : float
        Third rotation angle (according to `axes`).
    Examples
    --------
    >>> R0 = euler2mat(1, 2, 3, 'syxz')
    >>> al, be, ga = mat2euler(R0, 'syxz')
    >>> R1 = euler2mat(al, be, ga, 'syxz')
    >>> np.allclose(R0, R1)
    True
    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    M = np.array(mat, dtype=np.float64, copy=False)[:3, :3]
    if repetition:
        sy = np.sqrt(M[i, j]*M[i, j] + M[i, k]*M[i, k])
        if sy > _EPS4:
            ax = np.arctan2( M[i, j],  M[i, k])
            ay = np.arctan2( sy,       M[i, i])
            az = np.arctan2( M[j, i], -M[k, i])
        else:
            ax = np.arctan2(-M[j, k],  M[j, j])
            ay = np.arctan2( sy,       M[i, i])
            az = 0.0
    else:
        cy = np.sqrt(M[i, i]*M[i, i] + M[j, i]*M[j, i])
        if cy > _EPS4:
            ax = np.arctan2( M[k, j],  M[k, k])
            ay = np.arctan2(-M[k, i],  cy)
            az = np.arctan2( M[j, i],  M[i, i])
        else:
            ax = np.arctan2(-M[j, k],  M[j, j])
            ay = np.arctan2(-M[k, i],  cy)
            az = 0.0

    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax
    return ax/np.pi*180., ay/np.pi*180., az/np.pi*180.

def rot2mat(axis, angle, is_normalized=False):
    ''' Rotation matrix for rotation angle `angle` around `axis`
    ===============================================================================
    :::Important Note:::
    The rotation matrix is defined for rotation of a tensor in a fixed coordinate
    The output rotation matrix generated by this function is the inverse of the
    rotation matrix in Bond's book(p12-13).
    ===============================================================================
    Input Parameters:
    axis            - 3 element sequence, vector specifying axis for rotation.
    angle           - scalar, angle of rotation in degree.
    is_normalized   - bool, optional
       True if `axis` is already normalized (has norm of 1).  Default False
    -----
    output  -   mat : array shape (3,3), rotation matrix for specified rotation
    ===============================================================================
    Notes
    -----
    From: http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    '''
    x, y, z = axis
    if not is_normalized:
        n = np.sqrt(x*x + y*y + z*z)
        x = x/n
        y = y/n
        z = z/n
    angle  = np.pi*angle/180.
    c = np.cos(angle); s = np.sin(angle); C = 1-c
    xs = x*s;   ys = y*s;   zs = z*s
    xC = x*C;   yC = y*C;   zC = z*C
    xyC = x*yC; yzC = y*zC; zxC = z*xC
    return np.array([
            [ x*xC+c,   xyC-zs,   zxC+ys ],
            [ xyC+zs,   y*yC+c,   yzC-xs ],
            [ zxC-ys,   yzC+xs,   z*zC+c ]])

def mat2rot(mat, unit_thresh=1e-5):
    """Return axis, angle and point from (3, 3) matrix `mat`
    Parameters
    ----------
    mat : array-like shape (3, 3)
        Rotation matrix
    unit_thresh : float, optional
        Tolerable difference from 1 when testing for unit eigenvalues to
        confirm `mat` is a rotation matrix.
    Returns
    -------
    axis : array shape (3,)
       vector giving axis of rotation
    angle : scalar
       angle of rotation in radians.
    Examples
    --------
    >>> direc = np.random.random(3) - 0.5
    >>> angle = (np.random.random() - 0.5) * (2*math.pi)
    >>> R0 = axangle2mat(direc, angle)
    >>> direc, angle = mat2axangle(R0)
    >>> R1 = axangle2mat(direc, angle)
    >>> np.allclose(R0, R1)
    True
    Notes
    -----
    http://en.wikipedia.org/wiki/Rotation_matrix#Axis_of_a_rotation
    """
    M = np.asarray(mat, dtype=np.float)
    # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
    L, W = np.linalg.eig(M.T)
    i = np.where(np.abs(L - 1.0) < unit_thresh)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    direction = np.real(W[:, i[-1]]).squeeze()
    # rotation angle depending on direction
    cosa = (np.trace(M) - 1.0) / 2.0
    if abs(direction[2]) > 1e-8:
        sina = (M[1, 0] + (cosa-1.0)*direction[0]*direction[1]) / direction[2]
    elif abs(direction[1]) > 1e-8:
        sina = (M[0, 2] + (cosa-1.0)*direction[0]*direction[2]) / direction[1]
    else:
        sina = (M[2, 1] + (cosa-1.0)*direction[1]*direction[2]) / direction[0]
    angle = np.arctan2(sina, cosa)
    return direction, angle/np.pi*180.

def bondmat(axis, angle):
    """
    Compute Bond Matrix for rotation of Voigt matrix (eq. 8.9 in Bond, 1943; eq. 1.54 in Carcione, 2014)
    :::Important Note:::
    The rotation matrix used for Bond matrix was originally defined for rotation of the coordinate system
    ,which is in the opposite direction for rotation of a tensor in a fixed coordinate.
    We define the rotation as the rotation for the tensor itself in a fixed coordinate,
    Therefore,
    M   = bondmat(axis, angle) should be equal to the Bond Matrix for an angle of opposite sign
    =========================================================================================================
    Input Parameters:
    axis            - 3 element sequence, vector specifying axis for rotation.
    angle           - scalar, angle of rotation in degree.
    -----
    output          - array shape (3,3), Bond matrix for rotation of Voigt matrix
    =========================================================================================================
    """
    g       = rot2mat(axis = axis, angle = angle)
    M       = np.array([[g[0,0]**2, g[0,1]**2, g[0,2]**2, 2.*g[0,1]*g[0,2], 2.*g[0,2]*g[0,0], 2.*g[0,0]*g[0,1]],
                        [g[1,0]**2, g[1,1]**2, g[1,2]**2, 2.*g[1,1]*g[1,2], 2.*g[1,2]*g[1,0], 2.*g[1,0]*g[1,1]],
                        [g[2,0]**2, g[2,1]**2, g[2,2]**2, 2.*g[2,1]*g[2,2], 2.*g[2,2]*g[2,0], 2.*g[2,0]*g[2,1]],
        [g[1,0]*g[2,0], g[1,1]*g[2,1], g[1,2]*g[2,2], g[1,1]*g[2,2]+g[1,2]*g[2,1], g[1,0]*g[2,2]+g[1,2]*g[2,0], g[1,1]*g[2,0]+g[1,0]*g[2,1]],
        [g[2,0]*g[0,0], g[2,1]*g[0,1], g[2,2]*g[0,2], g[0,1]*g[2,2]+g[0,2]*g[2,1], g[0,2]*g[2,0]+g[0,0]*g[2,2], g[0,0]*g[2,1]+g[0,1]*g[2,0]],
        [g[0,0]*g[1,0], g[0,1]*g[1,1], g[0,2]*g[1,2], g[0,1]*g[1,2]+g[0,2]*g[1,1], g[0,2]*g[1,0]+g[0,0]*g[1,2], g[0,0]*g[1,1]+g[0,1]*g[1,0]]
        ])
    return M

def bondmat2(g):
    """
    Compute Bond Matrix for rotation of Voigt matrix (eq. 8.9 in Bond, 1943; eq. 1.54 in Carcione, 2014)
    ================================================================================
    Input Parameters:
    g   - transformation matrix
    -----
    output  -   mat : array shape (3,3), Bond matrix for rotation of Voigt matrix
    ================================================================================
    """
    M       = np.array([[g[0,0]**2, g[0,1]**2, g[0,2]**2, 2.*g[0,1]*g[0,2], 2.*g[0,2]*g[0,0], 2.*g[0,0]*g[0,1]],
                        [g[1,0]**2, g[1,1]**2, g[1,2]**2, 2.*g[1,1]*g[1,2], 2.*g[1,2]*g[1,0], 2.*g[1,0]*g[1,1]],
                        [g[2,0]**2, g[2,1]**2, g[2,2]**2, 2.*g[2,1]*g[2,2], 2.*g[2,2]*g[2,0], 2.*g[2,0]*g[2,1]],
        [g[1,0]*g[2,0], g[1,1]*g[2,1], g[1,2]*g[2,2], g[1,1]*g[2,2]+g[1,2]*g[2,1], g[1,0]*g[2,2]+g[1,2]*g[2,0], g[1,1]*g[2,0]+g[1,0]*g[2,1]],
        [g[2,0]*g[0,0], g[2,1]*g[0,1], g[2,2]*g[0,2], g[0,1]*g[2,2]+g[0,2]*g[2,1], g[0,2]*g[2,0]+g[0,0]*g[2,2], g[0,0]*g[2,1]+g[0,1]*g[2,0]],
        [g[0,0]*g[1,0], g[0,1]*g[1,1], g[0,2]*g[1,2], g[0,1]*g[1,2]+g[0,2]*g[1,1], g[0,2]*g[1,0]+g[0,0]*g[1,2], g[0,0]*g[1,1]+g[0,1]*g[1,0]]
        ])
    return M

def get_vel(bulk, shear, rho):
    """
    Compute primary and secondary seismic wave velocities for an isotropic material
    ================================================================================
    Input Parameters:
    bulk    - bulk modulus (GPa)
    shear   - shear modulus (Gpa)
    rho     - density (kg/m^3)
    -----
    output  - primary, secondary velocities (km/s)
    ================================================================================
    """
    primary     = np.sqrt(1000.0*(bulk + 4.0*shear/3)/rho)
    secondary   = np.sqrt(1000.0*shear/rho)
    return primary, secondary

def vec2rotmat(vector1, vector2):
    """
    Return a rotation matrix that rotates vector2 towards vector1.
    """
    vector1 = np.array(vector1)/norm(vector1)
    vector2 = np.array(vector2)/norm(vector2)
    rotvec  = np.cross(vector2, vector1)

    sin_angle = norm(rotvec)
    cos_angle = np.sqrt(1.0 - sin_angle*sin_angle)
    if sin_angle > 1e-10:
        dir_vec = rotvec/sin_angle
    else:
        return idmat

    ddt = np.outer(dir_vec, dir_vec)
    skew = np.array([[        0.0, -dir_vec[2],  dir_vec[1]],
                     [ dir_vec[2],         0.0, -dir_vec[0]],
                     [-dir_vec[1],  dir_vec[0],        0.0]])

    mtx = ddt + cos_angle * (idmat - ddt) - sin_angle * skew
    return mtx

def cofactor(m):
    """
    Return the cofactor matrix of a 3x3 matrix.
    """
    cof = np.empty((3, 3))

    cof[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1]
    cof[0][1] = m[1][2]*m[2][0] - m[1][0]*m[2][2]
    cof[0][2] = m[1][0]*m[2][1] - m[1][1]*m[2][0]

    cof[1][0] = m[0][2]*m[2][1] - m[0][1]*m[2][2]
    cof[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0]
    cof[1][2] = m[0][1]*m[2][0] - m[0][0]*m[2][1]

    cof[2][0] = m[0][1]*m[1][2] - m[0][2]*m[1][1]
    cof[2][1] = m[0][2]*m[1][0] - m[0][0]*m[1][2]
    cof[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0]
    
    cof2 = np.linalg.inv(m).T * np.linalg.det(m)
    if not np.allclose( cof, cof2):
        raise ValueError('Inconsistent cofactor matrix!')
    return cof

def localcartesian2spherical(xArr, yArr, zArr, thetaArr, phiArr, r=100):
    Ntheta, Nphi = thetaArr.shape
    uArr=np.zeros(thetaArr.shape)
    vArr=np.zeros(thetaArr.shape)
    rArr=np.zeros(thetaArr.shape)
    for itheta in xrange(Ntheta):
        for iphi in xrange(Nphi):
            theta   = thetaArr[itheta, iphi]/180.*np.pi
            phi     = phiArr[itheta, iphi]/180.*np.pi
            x       = xArr[itheta, iphi]
            y       = yArr[itheta, iphi]
            z       = zArr[itheta, iphi]
            if np.sin(theta) ==0: continue
            m=np.array([   [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(phi)],
                            [-np.sin(phi)/r/np.sin(theta), np.cos(phi)/r/np.sin(theta), 0],
                            [np.cos(phi)*np.cos(theta)/r, np.sin(phi)*np.cos(theta)/r, -np.sin(theta)/r]])
            out = np.dot(m, np.array([x,y,z]))
            rArr[itheta, iphi]=out[0]
            dphi = out[1]
            dtheta = out[2]
            l=np.sqrt(dphi**2+dtheta**2)
            
            uArr[itheta, iphi]=dphi/l
            vArr[itheta, iphi]=dtheta/l
    return rArr, uArr, vArr
            

class elasticTensor(object):
    """
    An object to manipulate elastic tensor in 3D coordinate
    ===========================================================================
    Cijkl   - 4th order elastic tensor (3*3*3*3, GPa)
    Cvoigt  - Voigt matrix (6*6, GPa)
    eCvoigt - error in Voigt matrix
    rho     - density (kg/m^3)
    compl   - element is compliance or not
    info    - auxliary information
    ===========================================================================
    """
    def __init__(self, Cvoigt = np.zeros([6,6]), eCvoigt = None, compl = False):
        self.Cijkl  = np.zeros([3,3,3,3])
        self.Cvoigt = Cvoigt
        self.eCvoigt= eCvoigt
        self.rho    = np.nan
        self.compl  = compl
        self.info   = ''
        if not np.allclose(Cvoigt, np.zeros([6,6])): self.Voigt2Cijkl()
        return
    
    def __str__(self):
        self.small2zero()
        outstr=self.info
        outstr=outstr+'\n------\nVoigt matrix (Gpa):'
        outstr=outstr+'\n'+self.Cvoigt.__str__()
        outstr=outstr+'\n------\ndensity = %g' %self.rho + ' km/m^3'
        return outstr
    
    def __repr__(self): return self.__str__()
    
    def small2zero(self, resetCijkl=True):
        self.Cvoigt[np.abs(self.Cvoigt)<self.Cvoigt.max()*1e-6]=0.
        if resetCijkl: self.Voigt2Cijkl()
        return
    
    def copy(self): return copy.deepcopy(self)
    
    def Cijkl2Voigt(self):
        """
        Convert 4th order elastic tensor to Voigt notation
        Use the optional argument "compl" for the elastic compliance (not 
        stiffness) tensor to deal with the multiplication 
        of elements needed to keep the Voigt and full 
        notation consistant.
        """
        t2m = np.array([[0,1,2,1,2,0],[0,1,2,2,0,1]])
        for i in xrange(6):
            for j in xrange(6):
                self.Cvoigt[i,j] = self.Cijkl[t2m[0,i],t2m[1,i],t2m[0,j],t2m[1,j]]
        if self.compl:
            self.Cvoigt = self.Cvoigt * np.array([  [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                    [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                    [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                    [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                                    [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                                    [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])
        return
    
    def Voigt2Cijkl(self):
        """
        Convert Voigt matrix to 4th order elastic tensor 
        Use the optional argument "compl" for the elastic compliance (not 
        stiffness) tensor to deal with the multiplication 
        of elements needed to keep the Voigt and full 
        notation consistant.
        """
        m2t = np.array([[0,5,4],[5,1,3],[4,3,2]])
        if self.compl:
            Cvoigt = self.Cvoigt / np.array([   [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                [1.0, 1.0, 1.0, 2.0, 2.0, 2.0],
                                                [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                                [2.0, 2.0, 2.0, 4.0, 4.0, 4.0],
                                                [2.0, 2.0, 2.0, 4.0, 4.0, 4.0]])
        else: Cvoigt = self.Cvoigt
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        self.Cijkl[i,j,k,l] = Cvoigt[m2t[i,j],m2t[k,l]]
        return
    
    def is_isotropic(self, tol=1e-5, verbose=True):
        """
        Check whether the elastic tensor is isotropic or not.
        """
        C = self.Cvoigt
        isISO = (abs(C[0,0]-C[1,1]) < tol) and (abs(C[0,0]-C[2,2]) < tol) and \
            (abs(C[0,1]-C[0,2]) < tol) and (abs(C[0,1]-C[1,2]) < tol) and \
            (abs(C[3,3]-C[4,4]) < tol) and (abs(C[3,3]-C[5,5]) < tol) and \
            (abs(C[0,3]) < tol) and (abs(C[0,4]) < tol) and (abs(C[0,5]) < tol) and \
            (abs(C[1,3]) < tol) and (abs(C[1,4]) < tol) and (abs(C[1,5]) < tol) and \
            (abs(C[2,3]) < tol) and (abs(C[2,4]) < tol) and (abs(C[2,5]) < tol) and \
            (abs(C[3,4]) < tol) and (abs(C[3,5]) < tol) and (abs(C[4,5]) < tol) and \
            (((C[0,0]-C[0,1])/2.0)-C[3,3] < tol)
        if verbose:
            if isISO: print 'Elastic tensor is isotropic !'
            else: print 'Elastic tensor is anisotropic !'
        return isISO
    
    def check_stability(self, verbose=True):
        """
        Check that the elastic constants matrix is positive definite 
        That is, check that the structure is stable to small strains. This
        is done by finding the eigenvalues of the Voigt elastic stiffness matrix
        by diagonalization and checking that they are all positive.
        See Born & Huang, "Dynamical Theory of Crystal Lattices" (1954) page 141.
        """
        (eigenvalues, eigenvectors) = np.linalg.eig(self.Cvoigt)
        if not (np.amin(eigenvalues) > 0.0):
            print 'Eigenvalues:', eigenvalues
            raise ValueError('Elastic tensor is not stable to small strains (Voigt matrix is not positive definite) !')
        if verbose: print 'Stability checked! Eigenvalues:', eigenvalues
        return 
    
    def check_symmetry(self, tol=1e-5, verbose=True):
        """
        Check that the elastic constants matrix is symmetric
        """
        Ctemp       = self.Cvoigt.copy()
        dC          = Ctemp.T - Ctemp
        maxdC       = (np.abs(dC)).max()
        if maxdC > tol: raise ValueError('Elastic tensor is not symmetric !')
        if verbose: print 'Symmetry checked! Maximum element in dC =', maxdC
        return
    
    ###############################################################
    # Methods for specifying elastic parameters
    ###############################################################
    
    def set_iso(self, vp, vs, rho, resetCijkl=True):
        """
        Set isotropic parameters given P and S wave velocity
        ============================================================================
        Input Parameters:
        vp, vs  - P/S wave velocity (km/s)
        rho     - Density (kg/m3) 
        resetCijkl  - reset 4th order tensor or not
        ============================================================================
        """
        vp  = vp*1e3
        vs  = vs*1e3
        self.Cvoigt[2,2] = vp*vp 
        self.Cvoigt[5,5] = vs*vs 
   
        self.Cvoigt[0,0] = self.Cvoigt[2,2]; self.Cvoigt[1,1] = self.Cvoigt[2,2] 
        self.Cvoigt[4,4] = self.Cvoigt[5,5]; self.Cvoigt[3,3] = self.Cvoigt[5,5]
        self.Cvoigt[0,1] = self.Cvoigt[2,2]-2.*self.Cvoigt[3,3]
        self.Cvoigt[0,2] = self.Cvoigt[0,1]; self.Cvoigt[1,2] = self.Cvoigt[0,1] 
        #  convert to GPa
        self.Cvoigt = self.Cvoigt*rho/1e9
        # make symmetrical
        for i in xrange(6):
             for j in xrange(6):
                 self.Cvoigt[j,i] = self.Cvoigt[i,j]
        self.rho    = rho
        if resetCijkl: self.Voigt2Cijkl()
        return
        
    def set_love(self, A, C, L, N, F, resetCijkl=True, mtype='VTI'):
        """
        Set Love parameters for a VTI media
        ============================================================================
        Input Parameters:
        A,C,L,N,F   - Love parameters (GPa)
        resetCijkl  - reset 4th order tensor or not
        ============================================================================
        """
        if mtype=='VTI':
            self.Cvoigt[:] =  np.array([[A, A-2.*N, F, 0., 0., 0.],
                                        [A-2.*N, A, F, 0., 0., 0.],
                                        [F, F, C, 0., 0., 0.],
                                        [0., 0., 0., L, 0., 0.],
                                        [0., 0., 0., 0., L, 0.],
                                        [0., 0., 0., 0., 0., N]])
            self.info='Love VTI'
        elif mtype=='HTI':
            self.Cvoigt[:] =  np.array([[C, F, F, 0., 0., 0.],
                                        [F, A, A-2*N, 0., 0., 0.],
                                        [F, A-2*N, A, 0., 0., 0.],
                                        [0., 0., 0., N, 0., 0.],
                                        [0., 0., 0., 0., L, 0.],
                                        [0., 0., 0., 0., 0., L]])
            self.info='Love HTI'
        if resetCijkl: self.Voigt2Cijkl()
        return
    
    def set_radial(self, vsv=3.57, vsh=3.74, vpv=6.14, vph=6.52, eta=0.87, rho=2790, resetCijkl=True):
        """
        Set Love parameters for a VTI media, given radial anisotropic S/P wave velocity
        Default values are the values from Point A, group 1 in Table 1. of Xie et al., 2015
        ====================================================================================
        Input Parameters:
        vsv         - SV wave velocity (km/s)
        vsh         - SH wave velocity (km/s)
        vpv         - PV wave velocity (km/s)
        vph         - PH wave velocity (km/s)
        eta         - F/(A-2*L)
        rho         - density (kg/m^3)
        resetCijkl  - reset 4th order tensor or not
        ====================================================================================
        Reference:
        Xie, J., Ritzwoller, M.H., Brownlee, S.J. and Hacker, B.R., 2015.
            Inferring the oriented elastic tensor from surface wave observations: preliminary application across the western United States.
            Geophysical Journal International, 201(2), pp.996-1021.
        """
        self.rho    = rho
        A           = rho*(vph**2)/1000.
        C           = rho*(vpv**2)/1000.
        N           = rho*(vsh**2)/1000.
        L           = rho*(vsv**2)/1000.
        F           = eta*(A-2*L)
        self.set_love(A=A, C=C, L=L, N=N, F=F, resetCijkl=resetCijkl)
        return
    
    def set_radial2(self, vp, vs, rho, xi, phi, eta, resetCijkl=True):
        """
        Output the elastic tensor given a set of radial anisotropy parameters
        as used typically in global seismology.  Average velocities are given by:
            15*rho*<Vp>^2 = 3*C + (8 + 4*eta)*A + 8*(1 - eta)*L
            15*rho*<Vs>^2 =   C + (1 - 2*eta)*A + (6 + 4*eta)*L + 5*N
        ============================================================================
        Input Parameters:
        vp  - Voigt average P wave velocity
        vs  - Voigt average shear wave velocity
        rho - Density
        xi  - (Vsh^2/Vsv^2) of horizontal waves
        phi - (Vpv^2/Vph^2)
        eta - C13/(C11 - 2*C44)
        ============================================================================
        """
        vp=vp*1e3
        vs=vs*1e3
        L = 15.*rho*((3.*phi + 8. + 4.*eta)*vs**2 - (phi + 1. - 2.*eta)*vp**2) \
                /((6. + 4.*eta + 5.*xi)*(3.*phi + 8. + 4.*eta) 
                 - 8.*(phi + 1. - 2.*eta)*(1. - eta)) 

        A = (15.*rho*vp**2 - 8.*(1. - eta)*L) / (3.*phi + 8. + 4.*eta) 
     
        F = eta*(A - 2.*L) 
        C = phi*A 
        N = xi*L 
        C12 = A - 2.*N
        self.Cvoigt[:]  =  np.array([[A, C12, F, 0., 0., 0.],
                                    [C12, A, F, 0., 0., 0.],
                                    [F, F, C, 0., 0., 0.],
                                    [0., 0., 0., L, 0., 0.],
                                    [0., 0., 0., 0., L, 0.],
                                    [0., 0., 0., 0., 0., N]])
        self.Cvoigt     = self.Cvoigt/1e9
        self.info       = 'radial VTI'
        if resetCijkl: self.Voigt2Cijkl()
        return
    
    def set_thomsen(self, vp, vs, eps, gamma, delta, rho, resetCijkl=True):
        """
        Set Thomsen parameters for a VTI media
        ============================================================================
        Input Parameters:
        vp,vs                   - (km/s)
        eps,gamma,delta         - Thomsen parameters, dimensionless
                                    eps     : P wave anisotropy
                                    gamma   : SH wave anisotropy
                                    delta   : elliptical anisotropy: eps = delta
        rho                     - density (kg/m^2)
        resetCijkl              - reset 4th order tensor or not
        ============================================================================
        """
        #  convert to m/s
        vp      = vp*1e3
        vs      = vs*1e3
        C       = np.zeros([6,6])
        C[2,2]  = vp*vp  # Eq 9a in Thomsen paper.
        C[3,3]  = vs*vs  # 9b
        C[5,5]  = C[3,3]*(2.0*gamma +1.0) # 8b
        C[0,0]  = C[2,2]*(2.0*eps +1.0) # 8a
        btm     = 2.0*C[3,3]
        term    = C[2,2] - C[3,3]
        ctm     = C[3,3]*C[3,3] - (2.0*delta*C[2,2]*term + term*term) 
        dsrmt   = (btm*btm - 4.0*ctm) 
        if dsrmt < 0: raise ValueError('S-velocity too high or delta too negative for Thomsen routine.')
        C[0,2]  = -btm/2.0 + np.sqrt(dsrmt)/2.0 # Eq 17
        C[0,1]  = C[0,0] - 2.0*C[5,5] 
        C[1,2]  = C[0,2] 
        C[4,4]  = C[3,3] 
        C[1,1]  = C[0,0] 
        # make symmetrical
        for i in xrange(6):
            for j in xrange(6):
                C[j,i] = C[i,j]
        #  convert to GPa
        C           = C*rho/1e9
        # output data
        self.Cvoigt = C
        self.rho    = rho
        self.info   = 'Thomsen VTI'
        if resetCijkl: self.Voigt2Cijkl()
        return
    
         
    def elastic_DB(self, mtype, resetCijkl=True):
        """
        Get elastic constant from predefined data base
        ============================================================================
        Input Parameters:
        mtype       - mineral type
        resetCijkl  - reset 4th order tensor or not
        ============================================================================
        """
        mtype   = mtype.lower()
        if mtype == 'olivine' or mtype == 'ol':
            self.info = 'Single crystal olivine (Abramson et al, JGR, 1997; doi:10.1029/97JB00682)'  
            self.Cvoigt[:] =  np.array([[320.5, 68.1, 71.6, 0., 0., 0.],
                                        [68.1, 196.5, 76.8, 0., 0., 0.],
                                        [71.6, 76.8, 233.5, 0., 0., 0.],
                                        [0., 0., 0., 64.0, 0., 0.],
                                        [0., 0., 0., 0., 77.0, 0.],
                                        [0., 0., 0., 0., 0., 78.7]])
            self.rho=3355.
        elif mtype == 'fayalite' or mtype == 'fa':
            self.info = 'Isentropic elastic constants for single cystal (Fe0.94,Mn0.06)2SiO4 fayalite (Speziale et al, JGR, 2004; doi:10.1029/2004JB003162)' 
            self.Cvoigt[:] =  np.array([[270., 103.,  97.2, 0., 0., 0.],
                                        [103., 171.1, 93.5, 0., 0., 0.],
                                        [97.2, 93.5, 234.1, 0., 0., 0.],
                                        [0., 0., 0., 33.4, 0., 0.],
                                        [0., 0., 0., 0., 48.7, 0.],
                                        [0., 0., 0., 0., 0., 59.6]])
            self.rho=4339.1
        elif mtype == 'albite' or mtype == 'alb':
            self.info = 'Single crystal albite (Brown et al, PCM, 2006; doi:10.1007/s00269-006-0074-1)' 
            self.Cvoigt[:] =  np.array([[69.9,  34.0,  30.8,   5.1, -2.4, -0.9],
                                        [34.0, 183.5,   5.5,  -3.9, -7.7, -5.8],
                                        [30.8,   5.5, 179.5,  -8.7,  7.1, -9.8],
                                        [ 5.1,  -3.9,  -8.7,  24.9, -2.4, -7.2],
                                        [-2.4,  -7.7,   7.1,  -2.4, 26.8,  0.5],
                                        [-0.9,  -5.8,  -9.8,  -7.2,  0.5,  33.5]])
            self.rho=2623.
        elif mtype == 'anorthite' or mtype == 'an96':
            # Brown, Angel and Ross "Elasticity of Plagioclase Feldspars" 
            self.info = 'Single crystal anorthite 96 (Brown et al, JGR, 2016)'
            self.Cvoigt[:] =  np.array([[132.2, 64.0,  55.3,   9.5,  5.1, -10.8],
                                        [ 64.0,200.2,  31.9,   7.5,  3.4,  -7.2],
                                        [ 55.3, 31.9, 163.9,   6.6,  0.5,   1.6],
                                        [  9.5,  7.5,   6.6,  24.6,  3.0,  -2.2],
                                        [  5.1,  3.4,   0.5,   3.0, 36.6,   5.2],
                                        [-10.8, -7.2,   1.6,  -2.2,  5.2,  36.0]])
            self.rho=2757.
        elif mtype == 'enstatite' or mtype == 'ens':
            self.info = 'Single crystal orthoenstatite (Weidner et al, PEPI 1978, 17:7-13)'
            self.Cvoigt[:] =  np.array([[225.0, 72.0, 54.0, 0., 0., 0.],
                                        [ 72.0,178.0, 53.0, 0., 0., 0.],
                                        [ 54.0, 53.0,214.0, 0., 0., 0.],
                                        [0., 0., 0., 78.0, 0., 0.],
                                        [0., 0., 0., 0., 76.0, 0.],
                                        [0., 0., 0., 0., 0., 82.0]])
            self.rho=3200.
        elif mtype == 'jadeite' or mtype == 'jd':
            self.info = 'Single crystal jadeite (Kandelin and Weidner, 1988, 50:251-260)'
            self.Cvoigt[:] =  np.array([[274.0, 94.0, 71.0, 0., 4., 0.],
                                        [ 94.0,253.0, 82.0, 0.,14., 0.],
                                        [ 71.0, 82.0,282.0, 0.,28., 0.],
                                        [   0.,   0.,   0.,88.0, 0., 13.0],
                                        [  4.0,  14.,  28.,  0.,65.,  0.],
                                        [   0.,   0.,   0., 13.0, 0., 94.0]])
            self.rho=3300.
        elif mtype == 'diopside' or mtype == 'di':
            self.info = 'Single crystal chrome-diopeside (Isaak and Ohno, PCM, 2003, 30:430-439)'
            self.Cvoigt[:] =  np.array([[228.1, 94.0, 71.0, 0., 7.9, 0.],
                                        [ 94.0,181.1, 82.0, 0., 5.9, 0.],
                                        [ 71.0, 82.0,245.4, 0.,39.7, 0.],
                                        [   0.,  0.0,  0.0,78.9, 0.,6.4],
                                        [  7.9,  5.9, 39.7,  0.,68.2,0.],
                                        [   0.,   0.,   0., 6.4, 0., 78.1]])
            self.rho=3400.
        elif mtype == 'halite' or mtype == 'nacl':
            self.info = 'Single crystal halite (NaCl, rock-salt).'
            self.Cvoigt[:] =  np.array([[ 49.5, 13.2, 13.2, 0.0, 0., 0.],
                                        [ 13.2, 49.5, 13.2, 0.0, 0., 0.],
                                        [ 13.2, 13.2, 49.5, 0.0, 0., 0.],
                                        [   0.,  0.0,  0.0,12.8, 0., 0.],
                                        [   0.,   0.,   0.,  0.,12.8,0.],
                                        [   0.,   0.,   0.,  0., 0., 12.8]])
            self.rho=2170.
        elif mtype == 'sylvite' or mtype == 'kcl':
            self.info = 'Single crystal sylvite (KCl).'
            self.Cvoigt[:] =  np.array([[ 40.1,  6.6,  6.6, 0.0, 0., 0.],
                                        [  6.6, 40.1,  6.6, 0.0, 0., 0.],
                                        [  6.6,  6.6, 40.1, 0.0, 0., 0.],
                                        [   0.,  0.0,  0.0, 6.4, 0., 0.],
                                        [   0.,   0.,   0.,  0.,6.4, 0.],
                                        [   0.,   0.,   0.,  0., 0., 6.4]])
            self.rho=1990.
        elif mtype == 'galena':
            self.info = 'Single crystal galena (Bhagavantam and Rao, Nature, 1951 168:42)'
            self.Cvoigt[:] =  np.array([[ 127., 29.8, 29.8, 0.0, 0., 0.],
                                        [ 29.8, 127., 29.8, 0.0, 0., 0.],
                                        [ 29.8, 29.8, 127., 0.0, 0., 0.],
                                        [   0.,  0.0,  0.0,24.8, 0., 0.],
                                        [   0.,   0.,   0.,  0.,24.8, 0.],
                                        [   0.,   0.,   0.,  0., 0., 24.8]])
            self.rho=7600.
        elif mtype == 'stishovite':
            self.info = 'Single crystal stishovite, SiO2 (Weidner et al., JGR, 1982, 87:4740-4746)'
            self.Cvoigt[:] =  np.array([[ 453., 211., 203., 0.0, 0., 0.],
                                        [ 211., 453., 203., 0.0, 0., 0.],
                                        [ 203., 203., 776., 0.0, 0., 0.],
                                        [   0.,  0.0,  0.0,252., 0., 0.],
                                        [   0.,   0.,   0.,  0.,252., 0.],
                                        [   0.,   0.,   0.,  0., 0., 252.]])
            self.rho=4290.
        elif mtype == 'fluorapatite'or mtype == 'apatite':
            self.info = 'Single crystal fluorapatite, Ca5F(PO4)3 (Sha et al., J. Appl. Phys., 1994, 75:7784; doi:10.1063/1.357030)'
            self.Cvoigt[:] =  np.array([[ 152., 49.99,  63.11, 0.0, 0., 0.],
                                        [49.99, 152.0,  63.11, 0.0, 0., 0.],
                                        [63.11, 63.11,  185.7, 0.0, 0., 0.],
                                        [   0.,  0.0,  0.0,  42.75, 0., 0.],
                                        [   0.,   0.,   0.,  0., 42.75, 0.],
                                        [   0.,   0.,   0.,  0.,    0., 51.005]])
            self.rho=3150.
        elif mtype == 'antigorite' or mtype == 'atg':
            # Note that X1||a, X2||b and X3||c* - not IRE convection.
            # and that these are adiabatic, not isothermal, constants, but
            # that's what "should" be used for wave velocites (c.f. Karato
            # deformation of Earth materials, sec. 4.3). Velocities quoted
            # in the reference (Table 3) use the corrected isotropic 
            # adiabatic moduli... 
            self.info = 'Adiabatic single crystal antigorite, (Bezacier et al., EPSL 2010, 289:198-208; doi:10.1016/j.epsl.2009.11.009)'
            self.Cvoigt[:] =  np.array([[ 208.1,  66.4,  16.00, 0.0,  5.5, 0.],
                                        [  66.4, 201.6,   4.90, 0.0, -3.1, 0.],
                                        [ 16.00,  4.90,  96.90, 0.0,  1.6, 0.],
                                        [    0.,    0.0,   0.0,16.9,  0.0, -12.10],
                                        [   5.5,  -3.10,   1.6,  0., 18.4, 0.],
                                        [   0.0,    0.0,   0.0,-12.10, 0., 65.50]])
            self.rho=2620.
        elif mtype == 'llm_mgsio3ppv':
            self.info = 'Adiabatic single crystal MgSiO3 post-perovskite under lowermost mantle conditions: from molecular dynamics, DFT & GGA and 2800 K and 127 GPa \
                        (Wookey et al. Nature, 2005, 438:1004-1007; doi:10.1038/nature04345 )'
            self.Cvoigt[:] =  np.array([[ 1139.,  357.,  311., 0.0,   0., 0.],
                                        [  357.,  842.,  466., 0.0,   0., 0.],
                                        [  311.,  466., 1137., 0.0,   0., 0.],
                                        [    0.,    0.,   0.0, 268.,  0., 0.],
                                        [    0.,    0.,   0.0,  0., 210., 0.],
                                        [    0.,    0.,   0.0,  0.,   0., 346.]])
            self.rho=5269.
        elif mtype == 'llm_mgsio3pv':
            self.info = 'Adiabatic single crystal MgSiO3 post-perovskite under lowermost mantle conditions: from molecular dynamics, DFT & GGA and 2800 K and 127 GPa \
                        (Wookey et al. Nature, 2005, 438:1004-1007; doi:10.1038/nature04345 )'
            self.Cvoigt[:] =  np.array([[  808.,  522.,  401., 0.0,   0., 0.],
                                        [  522., 1055.,  472., 0.0,   0., 0.],
                                        [  401.,  472.,  993., 0.0,   0., 0.],
                                        [    0.,    0.,   0.0, 328.,  0., 0.],
                                        [    0.,    0.,   0.0,  0., 263., 0.],
                                        [    0.,    0.,   0.0,  0.,   0., 262.]])
            self.rho=5191.
        elif mtype == 'ice':
            self.info = 'Adiabatic single crystal artifical water ice at -16C: from Gammon et al. (1983) Journal of Glaciology 29:433-460.'
            self.Cvoigt[:] =  np.array([[13.961,  7.153, 5.765, 0.0,   0., 0.],
                                        [ 7.153, 13.961, 5.765, 0.0,   0., 0.],
                                        [ 5.765,  5.765,15.013, 0.0,   0., 0.],
                                        [    0.,    0.,   0.0, 3.21,   0., 0.],
                                        [    0.,    0.,   0.0,   0., 3.21, 0.],
                                        [    0.,    0.,   0.0,   0.,   0., 3.404]])
            self.rho=919.10
        elif mtype == 'quartz' or mtype == 'qz':
            self.info = 'Premium cultured single crystal alpha-quartz at 22 C using resonance-ultrasound spectroscopy. \
                        From Heyliger et al. (2003) Journal of the Accoustic Society of America 114:644-650'
            self.Cvoigt[:] =  np.array([[ 87.17,   6.61, 12.02, -18.23,    0., 0.],
                                        [  6.61,  87.17, 12.02,  18.23,    0., 0.],
                                        [ 12.02,  12.02, 105.8,    0.0,    0., 0.],
                                        [-18.23,  18.23,   0.0,  58.27,    0., 0.],
                                        [    0.,     0.,   0.0,     0., 58.27, -18.23],
                                        [    0.,     0.,   0.0,     0.,-18.23,  40.28]])
            self.rho=2649.7
        elif mtype == 'periclase':
            self.info = 'Karki, B., Wentzcovitch, R., Gironcoli, S.D. & Baroni, S., 2000. Ab initiolattice dynamics of MgSiO3 '+\
                        'perovskite at high pressure, Phys. Rev. B,62(22), doi:10.1103/PhysRevB.62.14750'
            self.Cvoigt[:] =  np.array([[ 1154.,   265.5, 265.5, 0.,    0., 0.],
                                        [ 265.5,  1154., 265.5,  0.,    0., 0.],
                                        [ 265.5,  265.5, 1154.,  0.,    0., 0.],
                                        [   0.0,    0.0,   0.0, 198.0,  0., 0.],
                                        [    0.,     0.,   0.0,     0., 198.0, 0.],
                                        [    0.,     0.,   0.0,     0.,     0.,   198.0]])
            self.rho=5070.
        elif mtype == 'perovskite':
            self.info = 'Wentzcovitch, R., 2010. Thermodynamic properties and phase relations in mantle minerals investigated by'+\
                        'first principles quasiharmonic theory, Rev. Miner. Geochem., 71(1), 59?98.'
            self.Cvoigt[:] =  np.array([[ 860.,   535.5, 437.0, 0.,    0., 0.],
                                        [ 535.5,  1067.5,467.5,  0.,    0., 0.],
                                        [ 437.0,  467.5, 1053.,  0.,    0., 0.],
                                        [   0.0,    0.0,   0.0, 294.,  0., 0.],
                                        [    0.,     0.,   0.0,     0., 249.5, 0.],
                                        [    0.,     0.,   0.0,     0.,     0.,   284.5]])
            self.rho=5250.
        elif mtype == 'post-perovskite':
            self.info = 'Stackhouse, S. & Brodholt, J., 2007. The high-temperature elasticity of MgSiO3 post-perovskite, Geophys.'+\
                        'Monogr. Ser., 174, 99?114.'
            self.Cvoigt[:] =  np.array([[ 1220.,   474., 359., 0.,    0., 0.],
                                        [ 474.,  899., 493.,  0.,    0., 0.],
                                        [ 359.,  493., 1176.,  0.,    0., 0.],
                                        [   0.0,    0.0,   0.0, 273.0,  0., 0.],
                                        [    0.,     0.,   0.0,     0., 245.0, 0.],
                                        [    0.,     0.,   0.0,     0.,     0.,   376.0]])
            self.rho=5350.
        else: raise NameError('Unexpected name of mineral !')
        if resetCijkl: self.Voigt2Cijkl()
        return
    
    ##########################################################################
    # Methods for transformation of elastic tensor in fixed coordinate system
    ##########################################################################
    
    def rotT(self, axis, angle, resetCvoigt=True):
        """
        Rotate a 4th order elastic tensor with transformation matrix (rotation matrix),
        the coordinate system is fixed. 
        Note that the rotation is the inverse of rotation of a coordinate system.
        ==================================================================================
        Input Parameters:
        axis            - 3 element sequence, vector specifying axis for rotation.
        angle           - scalar, angle of rotation in degree.
        ==================================================================================
        ::: Note 1 :::
        Define L[i,j] = e'[i] * e[j] as transformation matrix of the coordinate,
        for axis = [0, 0, 1], rotation angle theta, we get:
            L = [cos(theta)  sin(theta)   0
                 -sin(theta) cos(theta)   0
                    0           0         1 ]
        And C'[i,j,k,l] = L[i,a]L[j,b]L[k,c]L[l,d]C[a,b,c,d]
        Note that we actually need to rotate the tensor, thus the angle has opposite sign,
        which is the case for the output from rot2mat.
        For more details, see Riley's book(p931, third edition).
        ::: Note 2 :::
        The definition of L[i,j] may be different in other books.
        Another definition is L[i,j] = e[i] * e'[j], in this case, we have:
        C'[i,j,k,l] = L[a,i]L[b,j]L[c,k]L[d,l]C[a,b,c,d]
        The rotation matrix from rot2mat should be changed if using this convention.
        ==================================================================================
        """
        g  = rot2mat(axis=axis, angle=angle)
        if use_opt_einsum:
            self.Cijkl[:]=contract('ia,jb,kc,ld,abcd->ijkl', g, g, g, g, self.Cijkl)
        else:
            self.Cijkl[:]=np.einsum('ia,jb,kc,ld,abcd->ijkl', g, g, g, g, self.Cijkl)
        if resetCvoigt: self.Cijkl2Voigt()
        return
    
    def _test_rotT(self, v=np.array([1.,0.,0.]), axis=[0,0,1], angle=45.):
        R=pysat.rot2mat(axis=axis, angle=angle)
        vprime=np.einsum('ia,a->i', R, v)
        print 'Testing rotation of a vector with fixed coordinate'
        print 'v =',v,' vprime = ',vprime
        return
    
    def rotB(self, axis, angle, resetCijkl=True):
        """
        Rotate Voigt matrix using Bond matrix (eq. 1.58 in Carcione, 2014)
        Note that the rotation is the inverse of rotation of a coordinate system,
        thus the rotation matrix used to construct Bond matrix is the inverse of the
        rotation matrix in Bond's book (p12-13)
        ============================================================================
        Input Parameters:
        axis            - 3 element sequence, vector specifying axis for rotation.
        angle           - scalar, angle of rotation in degree.
        resetCijkl      - reset 4th order tensor or not
        ============================================================================
        """
        M           = bondmat(axis=axis, angle=angle)
        self.Cvoigt = np.dot(M, self.Cvoigt)
        self.Cvoigt = np.dot(self.Cvoigt, M.T)
        if resetCijkl: self.Voigt2Cijkl()
        return
    
    def rotTB(self, axis, angle, verbose=True):
        """
        Rotate elastic tensor with both rotT and rotB, output error if incompatible
        ============================================================================
        Input Parameters:
        axis            - 3 element sequence, vector specifying axis for rotation.
        angle           - scalar, angle of rotation in degree.
        ============================================================================
        """
        et_temp = self.copy()
        self.rotT(axis=axis, angle=angle)
        et_temp.rotB(axis=axis, angle=angle)
        if not np.allclose(self.Cvoigt, et_temp.Cvoigt):
            raise ValueError('Inconsistent Rotation!')
        else:
            if verbose: print 'Consistent rotation!'
        return
    
    def rot_dip_strike(self, dip, strike, method='default'):
        """
        Rotate elastic tensor dip and strike angle, original tensor should be VTI
        Definition of geographical coordinate system:
        x   - North; y   - East; z  - depth
        ============================================================================
        Input Parameters:
        dip             - dip angle in degree (0<=dip<=90.)
        strike          - strike angle in degree (0<=strike<360.)
                            clockwise from North (x-axis)
        method          - 'default' : rotate with dip and strike in two steps
                          'euler'   : rotate with Euler angles
        ============================================================================
        """
        if dip >90. or dip < 0.: raise ValueError('Dip should be within [0., 90.]!')
        if method=='default':
            self.rotTB(axis=[1,0,0], angle=dip)
            self.rotTB(axis=[0,0,1], angle=strike)
        elif method == 'euler':
            g  = euler2mat(ai=dip, aj=0., ak=strike, axes='sxyz')
            if use_opt_einsum:
                self.Cijkl[:]=contract('ia,jb,kc,ld,abcd->ijkl', g, g, g, g, self.Cijkl)
            else:
                self.Cijkl[:]=np.einsum('ia,jb,kc,ld,abcd->ijkl', g, g, g, g, self.Cijkl)
            self.Cijkl2Voigt()
        return
    
    def rot_dip_strike2(self, dip, strike, verbose=True):
        et_temp = self.copy()
        self.rot_dip_strike(dip=dip, strike=strike)
        et_temp.rot_dip_strike(dip=dip, strike=strike, method='euler')
        if not np.allclose(self.Cvoigt, et_temp.Cvoigt):
            raise ValueError('Inconsistent dip/strike Rotation!')
        else:
            if verbose: print 'Consistent dip/strike Rotation!'
        return
    
    ##########################################################################
    # Advanced methods
    ##########################################################################
    
    def decompose_MN(self):
        """
        Decompose Voigt matrix into
        (1) ETI: an effective trasversely isotropic(azimuthally independent) part
        (2) AA: an azimuthally dependent part
        ============================================================================
        Output:
        etETI   - elastic tensor object for ETI part
        etAA    - elastic tensor object for AA part
        ============================================================================
        Reference:
        Montagner, J.P. and Nataf, H.C., 1986. A simple method for inverting the azimuthal anisotropy of surface waves.
            Journal of Geophysical Research: Solid Earth, 91(B1), pp.511-520.
        """
        etETI   = self.copy()
        etAA    = self.copy()
        Cvoigt  = self.Cvoigt.copy()
        A       = 3.*(Cvoigt[0,0] + Cvoigt[1,1])/8. + Cvoigt[0,1]/4. + Cvoigt[5,5]/2.
        C       = Cvoigt[2,2]
        N       = (Cvoigt[0,0] + Cvoigt[1,1])/8. - Cvoigt[0,1]/4. + Cvoigt[5,5]/2.
        L       = (Cvoigt[3,3] + Cvoigt[4,4]) / 2.
        F       = (Cvoigt[0,2] + Cvoigt[1,2])/2
        etETI.set_love(A = A, C = C, L=L, N=N, F=F, resetCijkl=True, mtype='VTI')
        self.etETI  = etETI
        etAA.Cvoigt = self.Cvoigt - etETI.Cvoigt
        etAA.Voigt2Cijkl()
        self.etAA   = etAA
        return 
    
    def set_error(self, eCvoigt=np.zeros([6,6])): self.eCvoigt= eCvoigt
        
    def invert(self, resetCijkl=True):
        """
        Given a square matrix and a square matrix of the errors
        on each element, return the inverse of the matrix and the 
        propogated errors on the inverse.
        We use numpy for the inversion and eq.10 of Lefebvre, 
        Keeler, Sobie and White ('Propagation of errors for 
        matrix inversion' Nuclear Instruments and Methods in 
        Physics Research A 451 pp.520-528; 2000) to calculate 
        the errors. The errors can be reported directly as the 
        errors on the inverse matrix but to do useful further 
        propogation we need to report the covar matrix too.
        This is calculated from eq.9 and we then extract the 
        diagonal elements to get the errors (rather than implementing
        eq.10 too).
        Tested with the matrix:
                0.700(7) 0.200(2)
                0.400(4) 0.600(6)               
        which gives back the inverse and squared errors reported
        in Table 1 of the above reference.
        This is coded up for an elastic constants matrix (Cij) and 
        its inverse (the elastic compliance matrix, Sij), but should
        work for any rank 2 square matrix.
        
        Need benchmark!
        """
        Cij  = self.Cvoigt
        eCij = self.eCvoigt
        if eCij is None: raise ValueError('Need to specify error matrix!')
        # Assuming we have a rank 2 square array
        # of the same size for input array. 
        if (np.ndim(Cij) != 2):
            raise ValueError, "Matrix must be rank 2"
        if (np.shape(Cij)[0] != np.shape(Cij)[1]):
            raise ValueError, "Matrix must be square"
        if (np.shape(Cij) != np.shape(eCij)):
            raise ValueError, "Matrix and error matrix must have same rank and shape"
        # Calculate the inverse using numpy
        Sij = np.linalg.inv(Cij)
        # Set up output arrays (init as zeros)
        eSij = np.zeros_like(eCij)
        array_size = eSij[0].size
        vcovSij = np.zeros((array_size,array_size,array_size,array_size),dtype=type(eSij))
        # Build covariance arrays (i.e COV(C^-1[a,b],S^-1[b,c] - a 4d array).
        # This is an implementation of eq.9 of Lefebvre et al.
        for a in xrange (array_size):
            for b in xrange (array_size):
                for c in xrange (array_size):
                    for d in xrange (array_size):
                        for i in xrange (array_size):
                            for j in xrange (array_size):
                                vcovSij[a,b,c,d] = vcovSij[a,b,c,d] + \
                                 ((Sij[a,i]*Sij[c,i]) * (eCij[i,j]**2) * (Sij[j,b]*Sij[j,d]))
        # Extract the "diagonal" terms, which are
        # the errors on the elements of the inverse
        # and could also be calculated using eq.10
        for a in xrange (array_size):
            for b in xrange (array_size):                
                eSij[a,b] = np.sqrt(vcovSij[a,b,a,b])
        self.Cvoigt     = Sij
        self.eCvoigt    = eSij
        self.vcovCvoigt = vcovSij
        self.compl      = not self.compl
        if resetCijkl: self.Voigt2Cijkl()
        return 
    
    def VRHavg(self):
        """
        Compute Voight-Reuss-Hill average of elastic constants tensor
        and propogated error given the 6*6 matrix of elastic constants
        and the 6*6 matrix of errors. The errors are optional. Assumes
        no co-varance between the errors on the elastic constants but 
        does include the covariance on the (calculated) compliance 
        matrix.
        Need benchmark!
        """
        if self.compl: raise ValueError('Elastic tensor is compliance!')
        Cij = self.Cvoigt
        eCij= self.eCvoigt
        # Need compliances too:
        if eCij is None:
            sij = np.linalg.inv(Cij)
        else:
            complTensor = self.copy()
            complTensor.invert()
            sij = complTensor.Cvoigt
            eSij= complTensor.eCvoigt
            covSij = complTensor.vcovCvoigt
        # These equations are valid for all crystal systems (only 9 of 
        # the 21 elastic constants ever have to be used, e.g. see Anderson 
        # theory of the Earth, pg. 122 or the introduction to Hill, 1952).
        voigtB = (1.0/9)*(Cij[0,0] + Cij[1,1] + Cij[2,2] ) \
               + (2.0/9)*(Cij[0,1] + Cij[0,2] + Cij[1,2])
        if eCij is not None:
            evB = np.sqrt( (1.0/81)*(eCij[0,0]**2 + eCij[1,1]**2 + eCij[2,2]**2) \
                      +(2.0/81)*(eCij[0,1]**2 + eCij[0,2]**2 + eCij[1,2]**2) )
        reussB = 1.0/((sij[0,0]+sij[1,1]+sij[2,2]) + 2*(sij[0,1]+sij[0,2]+sij[1,2]))
        if eCij is not None:
            # Note that COV(X+Z,Y) = COV(X,Y)+COV(Z,Y) and 
            # COV(SUM(Xi),SUM(Yj)) = SUM(SUM(COV(Xi,Yj)
            # c.f. http://mathworld.wolfram.com/Covariance.html
            erB = (np.sqrt(eSij[0,0]**2 + eSij[1,1]**2 + eSij[2,2]**2  \
                       + 4*eSij[0,1]**2 + 4*eSij[0,2]**2 + 4*eSij[1,2]**2  \
                       + 2*covSij[0,0,1,1] + 2*covSij[0,0,2,2] + 2*covSij[1,1,2,2] \
                       + 4*covSij[0,0,0,1] + 4*covSij[0,0,0,2] + 4*covSij[0,0,1,2] \
                       + 4*covSij[1,1,0,1] + 4*covSij[1,1,0,2] + 4*covSij[1,1,1,2] \
                       + 4*covSij[2,2,0,1] + 4*covSij[2,2,0,2] + 4*covSij[2,2,1,2] \
                       + 8*covSij[0,1,0,2] + 8*covSij[0,1,1,2] + 8*covSij[0,2,1,2] )) \
                * reussB**2
        voigtG = (1.0/15)*(Cij[0,0] + Cij[1,1] + Cij[2,2] - \
                           Cij[0,1] - Cij[0,2] - Cij[1,2]) + \
                 (1.0/5)*(Cij[3,3] + Cij[4,4] + Cij[5,5])
        if eCij is not None:
            evG = np.sqrt( (1.0/225)*(eCij[0,0]**2 + eCij[1,1]**2 + \
                                  eCij[2,2]**2 + eCij[0,1]**2 + \
                                  eCij[0,2]**2 + eCij[1,2]**2) + \
                        (1.0/25)*(eCij[3,3]**2 + eCij[4,4]**2 + eCij[5,5]**2) )
        reussG = 15.0/(4*(sij[0,0]+sij[1,1]+sij[2,2]) - \
                       4*(sij[0,1]+sij[0,2]+sij[1,2]) + 3*(sij[3,3]+sij[4,4]+sij[5,5]))
        if eCij is not None:
            erG = np.sqrt( \
                      16*(eSij[0,0]**2 + eSij[1,1]**2 + eSij[2,2]**2) \
                    + 16*(eSij[0,1]**2 + eSij[0,2]**2 + eSij[1,2]**2) \
                    +  9*(eSij[3,3]**2 + eSij[4,4]**2 + eSij[5,5]**2) \
                    + 32*covSij[0,0,1,1] + 32*covSij[0,0,2,2] + 32*covSij[1,1,2,2] \
                    + 32*covSij[0,0,0,1] + 32*covSij[0,0,0,2] + 32*covSij[0,0,1,2] \
                    + 32*covSij[1,1,0,1] + 32*covSij[1,1,0,2] + 32*covSij[1,1,1,2] \
                    + 32*covSij[2,2,0,1] + 32*covSij[2,2,0,2] + 32*covSij[2,2,1,2] \
                    + 32*covSij[0,1,0,2] + 32*covSij[0,1,1,2] + 32*covSij[0,2,1,2] \
                    + 24*covSij[0,0,3,3] + 24*covSij[0,0,4,4] + 24*covSij[0,0,5,5] \
                    + 24*covSij[1,1,3,3] + 24*covSij[1,1,4,4] + 24*covSij[1,1,5,5] \
                    + 24*covSij[2,2,3,3] + 24*covSij[2,2,4,4] + 24*covSij[2,2,5,5] \
                    + 24*covSij[0,1,3,3] + 24*covSij[0,1,4,4] + 24*covSij[0,1,5,5] \
                    + 24*covSij[0,2,3,3] + 24*covSij[0,2,4,4] + 24*covSij[0,2,5,5] \
                    + 24*covSij[1,2,3,3] + 24*covSij[1,2,4,4] + 24*covSij[1,2,5,5] \
                    + 18*covSij[3,3,4,4] + 18*covSij[3,3,5,5] + 18*covSij[4,4,5,5] \
                    ) * (reussG**2 / 15)
        if eCij is not None:
            return (voigtB, reussB, voigtG, reussG, ((voigtB+reussB)/2.0), ((voigtG+reussG)/2.0),
                   evB, erB, evG, erG, ((evB+erB)/2), ((evG+erG)/2))
        else:
            return (voigtB, reussB, voigtG, reussG, ((voigtB+reussB)/2.0), ((voigtG+reussG)/2.0),
                   None, None, None, None, None, None)

    def zenerAniso(self, eCvoigt=None):
        """
        Compute Zener anisotropy index, A, defined as
        2*C44/(C11-C12). This is unity for an isotropic crystal 
        and, for a cubic crystal C44 and 1/2(C11-C12) are shear 
        strains accross the (100) and (110) planes, respectivly.
        See Zener, Elasticity and Anelasticity of Metals, 1948
        or doi:10.1103/PhysRevLett.101.055504 (c.f. uAniso).
        Also returns the error on the anisotriopy index.
        Note that we don't check that the crystal is cubic!
        """
        Cij  = self.Cvoigt
        if eCvoigt is None: eCij = self.eCvoigt
        else: eCij = eCvoigt
        zA = (Cij[3,3]*2)/(Cij[0,0]-Cij[0,1])
        if eCij is None:
            return zA, None
        else:
            ezA = np.sqrt(((eCij[0,0]/Cij[0,0])**2 + (eCij[0,1]/Cij[0,1])**2) +\
               (2*(eCij[3,3]/Cij[3,3])**2)) * zA
            return (zA, ezA)

    def uAniso(self):
        """
        Returns the Universal elastic anisotropy index defined 
        by Ranganathan and Ostoja-Starzewski (PRL 101, 05504; 2008
        doi:10.1103/PhysRevLett.101.055504 ). Valid for all systems.
        """
        (voigtB, reussB, voigtG, reussG, hillB, hillG, 
                           evB, erB, evG, erG, ehB, ehG) = self.VRHavg()
        eCij= self.eCvoigt
        uA = (5*(voigtG/reussG))+(voigtB/reussB)-6
        if eCij is None:
            return uA, None
        else:
            euA = np.sqrt((np.sqrt((evG/voigtG)**2 + (erG/reussG)**2)*(voigtG/reussG))**2 + \
                      (np.sqrt((evB/voigtB)**2 + (erB/reussB)**2)*(voigtB/reussB))**2)
            return (uA, euA)

    def youngsmod(self, eCvoigt=np.zeros((6,6))):
        """
        Returns the Young's moduli, Poission ratio 
        and errors. Young's moduli is the ratio of tensile
        stress to tensile strain. Poission's ratios are ratio
        of tensile elongation and transverse contraction.
        """
        complTensor = self.copy()
        if complTensor.eCvoigt is None:  complTensor.set_error(eCvoigt=eCvoigt)
        complTensor.invert()
        sij     = complTensor.Cvoigt
        esij    = complTensor.eCvoigt
        covsij  = complTensor.vcovCvoigt
        
        youngX  = 1/sij[0,0]
        youngY  = 1/sij[1,1]
        youngZ  = 1/sij[2,2]
    
        eyoungX = (esij[0,0]/sij[0,0])*youngX
        eyoungY = (esij[1,1]/sij[1,1])*youngY
        eyoungZ = (esij[2,2]/sij[2,2])*youngZ
    
        poissonXY = -1*sij[0,1]*youngX
        poissonXZ = -1*sij[0,2]*youngX
        poissonYX = -1*sij[1,0]*youngY
        poissonYZ = -1*sij[1,2]*youngY
        poissonZX = -1*sij[2,0]*youngZ
        poissonZY = -1*sij[2,1]*youngZ
    
        epoissonXY = np.sqrt((esij[0,1]/sij[0,1])**2 + (esij[0,0]/sij[0,0])**2 - 
            2.0*((esij[0,1]*esij[0,0])/(sij[0,1]*sij[0,0]))*covsij[0,1,0,0])*poissonXY
        epoissonXZ = np.sqrt((esij[0,2]/sij[0,2])**2 + (esij[0,0]/sij[0,0])**2 - 
            2.0*((esij[0,2]*esij[0,0])/(sij[0,2]*sij[0,0]))*covsij[0,2,0,0])*poissonXZ
        epoissonYX = np.sqrt((esij[1,0]/sij[1,0])**2 + (esij[1,1]/sij[1,1])**2 - 
            2.0*((esij[1,0]*esij[1,1])/(sij[1,0]*sij[1,1]))*covsij[1,0,1,1])*poissonYX
        epoissonYZ = np.sqrt((esij[1,2]/sij[1,2])**2 + (esij[1,1]/sij[1,1])**2 - 
            2.0*((esij[1,2]*esij[1,1])/(sij[1,2]*sij[1,1]))*covsij[1,2,1,1])*poissonYZ
        epoissonZX = np.sqrt((esij[2,0]/sij[2,0])**2 + (esij[2,2]/sij[2,2])**2 - 
            2.0*((esij[2,0]*esij[2,2])/(sij[2,0]*sij[2,2]))*covsij[2,0,2,2])*poissonZX
        epoissonZY = np.sqrt((esij[2,1]/sij[2,1])**2 + (esij[2,2]/sij[2,2])**2 - 
            2.0*((esij[2,1]*esij[2,2])/(sij[2,1]*sij[2,2]))*covsij[2,1,2,2])*poissonZY
    
        return (youngX, youngY, youngZ, eyoungX, eyoungY, eyoungZ,
               poissonXY, poissonXZ, poissonYX, poissonYZ, poissonZX, poissonZY,
               epoissonXY, epoissonXZ, epoissonYX, epoissonYZ, epoissonZX, epoissonZY)

    def hessian_christoffelmat(self):
        """
        Return the hessian of the dynamical matrix.
        Due to the definition of the dynmat (q.C.q), this is independent of q.
        hessianmat[i][j][k][l] = d^2 M_kl / dx_i dx_j (note the indices).
        """
        hessianmat = np.empty((3, 3, 3, 3))
        for i in xrange(3):
            for j in xrange(3):
                for k in xrange(3):
                    for l in xrange(3):
                        hessianmat[i][j][k][l] = self.Cijkl[k][i][j][l] + self.Cijkl[k][j][i][l]
        return hessianmat
        
    def get_thomsen(self):
        """
        Get the Thomsen parameters.
        """
        C       = self.Cvoigt
        gamma   = (C[5,5] - C[3,3])/2./C[3,3]
        eps     = (C[0,0] - C[2,2])/2./C[2,2]
        delta   = ((C[0,2]+C[3,3])**2 - (C[2,2]-C[3,3])**2)/2./C[2,2]/(C[2,2]-C[3,3])
        return gamma, eps, delta
    
    def get_aziA(self, verbose=True):
        """
        Compute azimuthal anisotropy amplitude and corresponding fast axis angle
        """
        C       = self.Cvoigt
        try:
            etETI   = self.etETI
            L       = etETI.Cvoigt[4,4]
            if verbose: print   'Computing azimuthal anaisotropy WITH decomposed tensor!'
        except:
            L       = self.Cvoigt[4,4]
            if verbose: print   'Computing azimuthal anaisotropy WITHOUT decomposed tensor!'
        # # # aziA    = np.sqrt((C[4,4] - C[3,3])**2 +C[3,4]**2) / 2./L
        # # # phifa   = 1./2.*np.arctan(C[3,4]/ (C[4,4] - C[3,3])) / np.pi *180.
        aziA    = np.sqrt((C[4,4] - C[3,3])**2/4. +C[3,4]**2) / 2./L
        phifa   = 1./2.*np.arctan(2.*C[3,4]/ (C[4,4] - C[3,3])) / np.pi *180.
        return aziA, phifa
    
    def iso_vel(self, isotype='voigt'):
        VRHavg          = self.VRHavg()
        self.voigtB     = VRHavg[0]
        self.reussB     = VRHavg[1]
        self.voigtG     = VRHavg[2]
        self.reussG     = VRHavg[3]
        self.VRHB       = VRHavg[4]
        self.VRHG       = VRHavg[5]
        if isotype  =='voigt': iso_P, iso_S = get_vel(self.voigtB, self.voigtG, self.rho)
        elif isotype=='reuss': iso_P, iso_S = get_vel(self.reussB, self.reussG, self.rho)
        elif isotype=='VRH': iso_P, iso_S = get_vel(self.VRHB, self.VRHG, self.rho)
        return iso_P, iso_S
    
class Christoffel(object):
    """
    An object for solving (Kelvin-)Christoffel equation
    ===================================================================================
    Input Parameters:
    etensor     - elastic tensor object (type - elasticTensor)
    pv          - wave propagation direction vector
    -----------------------------------------------------------------------------------
    Output Parameters:
    kcmat       - Kelvin-Christoffel Matrix (3*3, unit: Gpa)
    grad_mat    - gradient of Kelvin-Christoffel Matrix (3*3*3, unit: Gpa)
    hessian_mat - Hessian matrix
    ----------------------
    --- Phase velocity ---
    ----------------------
    phvel       - phase velocities (3*1, unit: km/s)
    theta, phi  - polar angle, azimuth (scalar, unit: degree)
    eig_val     - eigenvalues of Kelvin-Christoffel equation (3*1, unit: GPa)
    eig_vec     - eigenvectors of Kelvin-Christoffel equation (polarization, 3*3)
    Note:   eig_vec[i, :] corresponds to eig_val[i]
    ----------------------
    --- Group velocity ---
    ----------------------
    grvel       - group velocities (3*1, unit: km/s)
    group_pv    - group velocity propagation vector (3*3)
    group_theta,- polar angle, azimuth of group velocity (3*1, unit: degree)
    group_phi  
    grad_eig_val- gradient of eigenvalues of Kelvin-Christoffel equation (3*3, unit: GPa)
    group_vec   - group velocity polorization vector (3*3)
    --------------
    --- Others ---
    --------------
    powflow_angle   - power flow angle (unit: degree)
                    (the angle between group/phase velocity propagation vector)
    cos_pf_angle    - cosine of power flow angle
    hessian_eig     - Hessian matrix of eigenvalues (3*3*3)
    ===================================================================================
    Modified from the Python code, christoffel by Jan Jaeken, added plotting function
    """

    def __init__(self, etensor, verbose=True, isotype='voigt'):
        # check input elastic tensor
        if not isinstance(etensor, elasticTensor): raise TypeError('Input object is not elasticTensor type!')
        etensor.check_stability(verbose=verbose)
        etensor.check_symmetry(verbose=verbose)
        etensor.is_isotropic(verbose=verbose)
        # Set values
        self.etensor    = etensor
        VRHavg          = etensor.VRHavg()
        self.voigtB     = VRHavg[0]
        self.reussB     = VRHavg[1]
        self.voigtG     = VRHavg[2]
        self.reussG     = VRHavg[3]
        self.VRHB       = VRHavg[4]
        self.VRHG       = VRHavg[5]
        if isotype  =='voigt': self.iso_P, self.iso_S = get_vel(self.voigtB, self.voigtG, self.etensor.rho)
        elif isotype=='reuss': self.iso_P, self.iso_S = get_vel(self.reussB, self.reussG, self.etensor.rho)
        elif isotype=='VRH': self.iso_P, self.iso_S = get_vel(self.VRHB, self.VRHG, self.etensor.rho)
        else: raise TypeError('Unexpected isotropic type !')
        self.hessian_mat = etensor.hessian_christoffelmat()
        self.clear_direction()
        return
    
    def copy(self): return copy.deepcopy(self)
            
    def clear_direction(self):
        """
        Clear all direction-dependent data
        """
        self.pv             = None
        self.theta          = None
        self.phi            = None
        self.kcmat          = None
        self.grad_mat       = None
        self.eig_val        = None
        self.eig_vec        = None
        self.phvel          = None
        self.grvel          = None
        self.group_vec      = None
        self.group_pv       = None
        self.group_theta    = None
        self.group_phi      = None
        self.grad_eig_val   = None
        self.powflow_angle  = None
        self.cos_pf_angle   = None
        self.hessian_eig    = None
        self.enhancement    = None
        return
    
    def set_direction_cartesian(self, pv, is_normalized=False):
        """
        Define a wave propagation vector in a Cartesian coordinate.
        It is always explicitly normalized to lie on the unit sphere.
        ============================================================================
        Input Parameters:
        pv              - wave propagation direction vector
                            list of 3 numbers or numpy array (e.g. [1, 0, 0])
        is_normalized   - whether the pv is unit vector or not
        ============================================================================
        """
        self.clear_direction()
        if not isinstance(pv, np.ndarray): pv = np.asarray(pv)
        x, y, z = pv
        if not is_normalized:
            n = np.sqrt(x*x + y*y + z*z)
            x = x/n
            y = y/n
            z = z/n
        if z == 1.0 or z == -1.0:
            if z > 0.0: self.theta = 0.0
            else: self.theta = np.pi
            self.phi = 0.0
        else:
            self.theta  = np.arccos(z)
            sin_theta   = np.sqrt(1 - z**2)
            cos_phi     = x/sin_theta
            self.phi    = np.arccos(cos_phi)
            if y < 0.0: self.phi = 2.0*np.pi - self.phi
        self.theta      = self.theta*180./np.pi
        self.phi        = self.phi*180./np.pi
        self.pv         = pv/n
        self.get_kc_mat()
        return
    
    def set_direction_spherical(self, theta, phi):
        """
        Define a wave vector in a spherical coordinate
        x = cos(phi) * sin(theta)
        y = sin(phi) * sin(theta)
        z = cos(theta)
        ==========================================================
        Input Parameters:
        theta   - polar angle (degree)
        phi     - azimuth (degree)
        ==========================================================
        """
        self.clear_direction()
        self.theta  = theta
        self.phi    = phi
        theta       = theta/180.*np.pi
        phi         = phi/180.*np.pi
        sin_theta   = np.sin(theta)
        cos_theta   = np.cos(theta)
        sin_phi     = np.sin(phi)
        cos_phi     = np.cos(phi)
        x           = cos_phi * sin_theta
        y           = sin_phi * sin_theta
        z           = cos_theta
        self.pv     = np.array([x, y, z])
        self.get_kc_mat()
        return
    
    def set_direction_random(self):
        """
        Generates a random wave vector direction.
        The distribution is uniform across the unit sphere.
        """
        self.clear_direction()
        cos_theta   = np.random.ranf()
        phi         = 2.0 * np.pi * np.random.ranf()
        sin_theta   = np.sqrt(1 - cos_theta**2)
        cos_phi     = np.cos(phi)
        sin_phi     = np.sin(phi)
        self.pv     = np.array([cos_phi*sin_theta, sin_phi*sin_theta, cos_theta])
        self.phi    = phi*180./np.pi
        self.theta  = np.arccos(cos_theta)*180./np.pi
        self.get_kc_mat()
        return
    
    def get_kc_mat(self):
        """
        Get the Kelvin-Christoffel Matrix (e.g. eq. 1.9 in Tsvankin, 2012)
        ===================================================================
        Output:
        self.kcmat  - Kelvin-Christoffel Matrix (3*3, unit: Gpa)
        ===================================================================
        """
        if use_opt_einsum:
            kcmat1  = contract('j,l, ijkl->ik', self.pv, self.pv, self.etensor.Cijkl)
        else:
            kcmat1  = np.einsum('j,l, ijkl->ik', self.pv, self.pv, self.etensor.Cijkl)
        ###
        # kcmat2  = np.dot(self.pv, np.dot(self.pv, self.etensor.Cijkl))
        # if not np.allclose(kcmat1, kcmat2): raise ValueError('Error Christoffel Matrix')
        ###
        self.kcmat=kcmat1
        return
    
    def get_phvel(self):
        """
        Determine eigenvalues, eigenvectors of the Kelvin-Christoffel matrix.
        Results are sorted from low to high, then store eigens and phase velocities.
        ========================================================================================
        Output:
        self.eig_val    - eigenvalues  (3*1, unit: Gpa)
        self.eig_vec    - eigenvectors (3*3, denote polirization direction of the wave)
        self.phvel      - phase velocity (3*1, unit: km/s)
        ========================================================================================
        """
        eig_val, eig_vec = np.linalg.eigh(self.kcmat)
        args    = np.argsort(eig_val)
        eig_val = eig_val[args]
        eig_vec = eig_vec.T[args]
        self.eig_val= eig_val
        self.eig_vec= eig_vec
        # eig_val has unit of Gpa
        self.phvel  = np.sign(eig_val)*np.sqrt(np.absolute(eig_val) * 1000. / self.etensor.rho )
        return
    
    def get_grad_mat(self):
        """
        Calculate the gradient of the Kelvin-Christoffel matrix.
        d/dx_k M_ij =  q_m * ( C_imkj + C_ikmj )
        gradmat[k][i][j] =  d/dx_k M_ij (note the indices)
        ========================================================================================
        Output:
        self.grad_mat  - gradient of the Kelvin-Christoffel Matrix (3*3*3, unit: Gpa)
        ========================================================================================
        """
        if use_opt_einsum:
            gradmat1 = contract('m, ikmj->kij', self.pv, self.etensor.Cijkl) + \
                    contract('m, imkj->kij', self.pv, self.etensor.Cijkl)
        else:
            gradmat1 = np.einsum('m, ikmj->kij', self.pv, self.etensor.Cijkl) + \
                    np.einsum('m, imkj->kij', self.pv, self.etensor.Cijkl)
        ###
        # # # pv  = self.pv
        # # # C   = self.etensor.Cijkl
        # # # gradmat2 = np.dot(pv, C + np.transpose(C, (0, 2, 1, 3)))
        # # # gradmat2 = np.transpose(gradmat2, (1, 0, 2))
        # # # if not np.allclose(gradmat1, gradmat2): raise ValueError('Error gradient Christoffel Matrix')
        ###
        self.grad_mat = gradmat1
        return
    
    def get_grvel(self):
        """
        Calculate group velocities as the gradient of the phase velocities.
        Power flow angles are also calculated and stored.
        ========================================================================================
        Output:
        self.grvel          - group velocities (3*1, unit: km/s)
        self.group_pv       - group velocity propagation vector 
        self.group_theta    - group velocity polar angles (3*1, unit: degree)
        self.group_phi      - group velocity azimuth (3*1, unit: degree)
        self.group_vec      - group velocity polorization vector (3*3)
        self.grad_eig_val   - gradient of eigenvalues (3*3, unit: GPa)
        self.cos_pf_angle   - cosine of power flow angle (3*1)
        self.powflow_angle  - power flow angle (3*1, unit: degree)
                                (the angle between group/phase velocity propagation vector)
        ========================================================================================
        """
        if self.phvel is None: self.get_phvel()
        if self.grad_mat is None: self.get_grad_mat()
        phase_vel   = self.phvel
        eig_vec     = self.eig_vec
        gradmat     = self.grad_mat * 1000. / self.etensor.rho 

        grad_eig        = np.empty((3, 3))
        group_vec       = np.empty((3, 3))
        self.grvel      = np.empty(3)
        self.group_pv   = np.empty((3, 3))
        self.group_theta= np.empty(3)
        self.group_phi  = np.empty(3)
        for pol in xrange(3):
            for cart in xrange(3):
                grad_eig[pol][cart] = \
                np.dot(eig_vec[pol], np.dot(gradmat[cart], eig_vec[pol]))
                # Eigenvalues are the square of the velocity
                # dv/dq = dv^2/dq / (2v)
                group_vec[pol][cart] = grad_eig[pol][cart] / (2*phase_vel[pol])
            self.grvel[pol] = np.linalg.norm(group_vec[pol], 2)
            self.group_pv[pol] = group_vec[pol] / self.grvel[pol]

            x = self.group_pv[pol][0]
            z = self.group_pv[pol][2]
            if z >= 1.0-1e-10 or z <= -1.0+1e-10:
                self.group_theta[pol] = 0.0
                self.group_phi[pol] = 0.0
            else:
                self.group_theta[pol] = np.arccos(z)
                sin_theta = np.sqrt(1 - z**2)
                if abs(x) > sin_theta:
                    self.group_phi[pol] = (1.0 - np.sign(x))*0.5*np.pi
                else:
                    self.group_phi[pol] = np.arccos(x/sin_theta)
                if self.group_pv[pol][1] < 0.0:
                    self.group_phi[pol] = 2*np.pi - self.group_phi[pol]
        # In case things go wrong, check if phase_vel == np.dot(group_vec, pv)
        self.grad_eig_val   = grad_eig * self.etensor.rho / 1000.
        self.group_vec      = group_vec 
        self.cos_pf_angle   = np.dot(self.group_pv, self.pv)
        self.powflow_angle  = np.arccos(np.around(self.cos_pf_angle, 10))/np.pi*180.
        self.group_theta    = self.group_theta*180./np.pi
        self.group_phi      = self.group_phi*180./np.pi
        if not np.allclose( self.phvel, np.dot(self.group_vec, self.pv)):
            raise ValueError('Inconsistent phase/group velocities!')
        return
    
    def get_hessian_eig(self):
        """
        Calculate the eigenvalues of Hessian matrix
        Hessian[n][i][j] = d^2 lambda_n / dx_i dx_j
        """
        dynmat  = self.kcmat
        eig_val = self.eig_val
        eig_vec = self.eig_vec
        gradmat = self.grad_mat
        hess_mat= self.hessian_mat
        hessian = np.zeros((3, 3, 3))
        idmat   = np.identity(3)
        for n in xrange(3):
            hessian[n]  += np.dot(np.dot(hess_mat, eig_vec[n]), eig_vec[n])
            pseudoinv   = np.linalg.pinv(eig_val[n]*idmat - dynmat, rcond=1e-10)
            deriv_vec   = np.dot(gradmat, eig_vec[n])
            hessian[n]  += 2.0 * np.dot(np.dot(deriv_vec, pseudoinv), deriv_vec.T)
            #Take deriv of eigenvec into account: 2 * (d/dx s_i) * pinv_ij * (d_dy s_j)
        self.hessian_eig= hessian
        return
    
    def get_enhancement(self, approx=False, num_steps=8, delta=1e-5):
        """
        Calculate the enhancement factor (p51, eq. 42 in Wolfe, 2005)
        ======================================================================================================
        Input Parameters:
        approx              - determine the enhancement factor approximately(numerically) or not
        self.enhancement    - enhancement factor
        num_steps           - number of sides for the surface polygon
        delta               - radius of the polygon
        ======================================================================================================
        When approx is True, the function determines the enhancement factors according to a numerical scheme.
        The surface areas of a set of triangles in phase and group space are
        calculated and divided. This is significantly slower and less accurate
        than the analytical approach, but will provide a physically relevant
        value when the enhancement factor is ill defined.

        The surface area is a polygon of n sides where n is num_steps.
        The radius of this polygon is determined by delta, which determines the
        change in theta and phi coordinates relative to the central position.
        
        Need benchmark!
        """
        if self.phvel is None: self.get_phvel()
        if self.grad_mat is None: self.get_grad_mat()
        if self.grvel is None: self.get_grvel()
        if self.hessian_eig is None: self.get_hessian_eig()
        if not approx:
            hessian     = self.hessian_eig *1000./self.etensor.rho
            phase_vel   = self.phvel
            group_vec   = self.group_vec
            group_abs   = self.grvel
            grad_group  = np.empty((3, 3, 3))
            enhance     = np.empty(3)
            for n in xrange(3):
                grad_group[n] = hessian[n] / group_abs[n]
                grad_group[n] -= np.outer(group_vec[n], np.dot(hessian[n], group_vec[n])) / (group_abs[n]**3)
                grad_group[n] /= 2.0*phase_vel[n] #grad lambda = 2 * v_p * v_g
    
                enhance[n] = 1.0 / np.linalg.norm( (np.dot(cofactor(grad_group[n]), self.pv)), 2)
            self.enhancement = enhance 
        else:
            tempChr     = self.copy()
            phase_grid  = np.empty((num_steps+1, 3))
            group_grid  = np.empty((num_steps+1, 3, 3))
    
            center_theta= self.theta/180.*np.pi
            center_phi  = self.phi/180.*np.pi
            phase_center= self.pv
    
            for i in xrange(num_steps):
                angle = i*2.0*np.pi/num_steps
                tempChr.set_direction_spherical( (center_theta + np.sin(angle)*delta)/np.pi*180., (center_phi + np.cos(angle)*delta)/np.pi*180.)
                phase_grid[i] = tempChr.pv
                if tempChr.group_pv is None: tempChr.get_grvel()
                group_grid[i] = tempChr.group_pv
    
            phase_grid[num_steps] = phase_grid[0]
            group_grid[num_steps] = group_grid[0]
            tempChr.set_direction_cartesian(phase_center)
            if tempChr.grvel is None: tempChr.get_grvel()
            group_center = tempChr.group_pv
    
            phase_area  = 0.0
            group_area  = np.zeros(3)
            tot_angle   = np.zeros(3)
            for i in xrange(num_steps):
                phase_area += np.linalg.norm ( (np.cross(phase_grid[i] - phase_center, phase_grid[i+1] - phase_center)), 2)
                for n in xrange(3):
                    group_area[n] += np.linalg.norm ( np.cross(group_grid[i][n] - group_center[n], group_grid[i+1][n] - group_center[n]), 2)
            self.enhancement = phase_area/group_area 
        return
    
    def find_nopowerflow(self, step_size=0.9, eig_id=2, max_iter=900):
        """
        Attempts to find the closest direction of extremal phase velocity,
        where group and phase directions align. A positive step_size should
        search for maxima, while negative step_size searches for minima.
        Due to the complicated nature of the ray surfaces of the quasi-shear
        modes (eig_id 0 and 1), there is no guarantee that this algorithm
        will converge or reliably find an extremal velocity.
        If a direction has been set already, the search will start from there
        and follow the general direction of power flow. Otherwise, the search
        will start from a randomly chosen point.
        """
        if self.pv is None:
            self.set_direction_random()
        phase_dir = self.pv
        if self.group_pv is None: self.get_grvel()
        group_dir = self.group_pv

        step_dir = group_dir[eig_id] - phase_dir
        if max_iter <= 0 or np.linalg.norm(step_dir, 2 ) < 1e-10:
            return
        else:
            self.set_direction_cartesian(phase_dir + step_size*step_dir)
            max_iter -= 1
            self.find_nopowerflow(step_size, eig_id, max_iter)
        return
    
    def circle(self, theta, dphi=1., group=False, verbose=False):
        """
        Scan a unit circle to produce data for azimuthal anisotropy variation figure
        ============================================================================================
        Input Parameters:
        theta   - polar angle
        dphi    - interval of azimuth
        group   - compute group velocities or nor
        --------------------------------------------------------------------------------------------
        Output:
        self.phvelArr       - phace velocity array (3*Nphi)
        self.eigvecArr      - eigenvector array (polarization of three phase modes, 3*3*Nphi)
        self.grvelArr       - group velocity array (3*Nphi)
        self.group_vecArr   - polarization of three group modes (3*3*Nphi)
        self.group_pvArr    - propagation vector of three group modes (3*3*Nphi)
        self.group_thetaArr - group velocity polar angle array (3*Nphi)
        self.group_phiArr   - group velocity azimuth array (3*Nphi)
        self.pf_angleArr    - power flow angle array (3*Nphi)
        ============================================================================================
        """
        Nphi    = int(360./dphi)
        phis    = np.linspace(0., 360., num=Nphi)
        # initialization
        phvelArr = np.zeros([3, Nphi])
        eigvecArr= np.zeros([3, 3, Nphi])
        if group:
            grvelArr        = np.zeros([3, Nphi])
            group_vecArr    = np.zeros([3, 3, Nphi])
            group_pvArr     = np.zeros([3, 3, Nphi])
            group_thetaArr  = np.zeros([3, Nphi])
            group_phiArr    = np.zeros([3, Nphi])
            pf_angleArr     = np.zeros([3, Nphi])
        print 'Solving Kelvin-Christoffel equation for a circle of propagation vector'
        dN= int(Nphi/10)
        i = 0; j = 0
        for iphi in xrange(Nphi):
            i+=1
            if i%dN ==0 and verbose:
                j+=10
                print 'Finished %'+str(j)
            phi         = phis[iphi]
            temp_kceq   = self.copy()
            temp_kceq.set_direction_spherical(theta=theta, phi=phi)
            temp_kceq.get_phvel()
            phvelArr[:, iphi]       = temp_kceq.phvel
            eigvecArr[:, :, iphi]   = temp_kceq.eig_vec
            if not group: continue
            temp_kceq.get_grvel()
            grvelArr[:, iphi]       = temp_kceq.grvel
            group_vecArr[:, :, iphi]= temp_kceq.group_vec
            group_pvArr[:, :, iphi] = temp_kceq.group_pv
            group_thetaArr[:, iphi] = temp_kceq.group_theta
            group_phiArr[:, iphi]   = temp_kceq.group_phi
            pf_angleArr[:, iphi]    = temp_kceq.powflow_angle
        print 'End solving Kelvin-Christoffel equation for a circle of propagation vector'
        self.phvelArr   = phvelArr
        self.eigvecArr  = eigvecArr
        self.phis       = phis
        if group:
            self.grvelArr        = grvelArr
            self.group_vecArr    = group_vecArr
            self.group_pvArr     = group_pvArr
            self.group_thetaArr  = group_thetaArr
            self.group_phiArr    = group_phiArr
            self.pf_angleArr     = pf_angleArr
        return
    
    def plot_circle(self, figsize=(15,15), polar=True, diff=True, rel=True, dtype='phase', showfig=True):
        if dtype is 'phase':
            R0  = self.phvelArr[0, :]
            R1  = self.phvelArr[1, :]
            R2  = self.phvelArr[2, :]
        elif dtype is 'group':
            R0  = self.grvelArr[0, :]
            R1  = self.grvelArr[1, :]
            R2  = self.grvelArr[2, :]
        else: raise NameError('Unexpected data type: ' + dtype)
        R0mean  = R0.mean()
        R1mean  = R1.mean()
        R2mean  = R2.mean()
        R0min  = R0.min()
        R1min  = R1.min()
        R2min  = R2.min()
        if diff:
            R0  = R0 - R0.min()
            R1  = R1 - R1.min()
            R2  = R2 - R2.min()
        # if rel:
        #     R0  = R0/R0mean*100.
        #     R1  = R1/R1mean*100.
        #     R2  = R2/R2mean*100.
        if rel:
            R0  = R0/R0min*100.
            R1  = R1/R1min*100.
            R2  = R2/R2min*100.
        phis= self.phis
        x0  = R0 * np.cos(np.pi/180.*phis)
        y0  = R0 * np.sin(np.pi/180.*phis)
        x1  = R1 * np.cos(np.pi/180.*phis)
        y1  = R1 * np.sin(np.pi/180.*phis)
        x2  = R2 * np.cos(np.pi/180.*phis)
        y2  = R2 * np.sin(np.pi/180.*phis)
        # plot results
        fig, ax = plt.subplots(figsize=figsize)
        if polar:
            ax = plt.subplot(111, projection='polar')
            plt.plot(phis/180*np.pi, R0, '-', lw=5, label='qS2 Wave')
            plt.plot(phis/180*np.pi, R1, '-', lw=5, label='qS1 Wave')
            plt.plot(phis/180*np.pi, R2, '-', lw=5, label='qP Wave')
        else:
            plt.plot(x0, y0, '-', lw=5, label='qS2 Wave')
            plt.plot(x1, y1, '-', lw=5, label='qS1 Wave')
            plt.plot(x2, y2, '-', lw=5, label='qP Wave')
            if rel:
                plt.ylabel('anisotropy(%)', fontsize=45)
                plt.xlabel('anisotropy(%)', fontsize=45)
            else:
                plt.ylabel('km/s', fontsize=45)
                plt.xlabel('km/s', fontsize=45)
            ax.tick_params(axis='y', labelsize=25)
            ax.tick_params(axis='x', labelsize=25)
            plt.axis('equal')
        plt.legend(numpoints=1, fontsize=20)
        plt.title(dtype+' velocity surface', fontsize=45)
        if rel: label = '%'
        else: label='km/s'
        print '=============================================================='
        print 'qP anisotropy: ', R2.max(),label
        print 'min vp =', R2min, 'km/s mean vp =',R2mean, 'km/s'
        print 'qS1 anisotropy: ', R1.max(),label
        print 'min vs1 =', R1min, 'km/s mean vs1 =',R1mean, 'km/s'
        print 'qS2 anisotropy: ', R0.max(),label
        print 'min vs2 =', R0min, 'km/s mean vs2 =',R0mean, 'km/s'
        print '=============================================================='
        if showfig: plt.show()
        return
    
    def sphere(self, dtheta=1., dphi=1., group=False, outfname=None):
        """
        Scan a unit sphere to produce data for pole figure
        ============================================================================================
        Input Parameters:
        dtheta  - interval of polar angle
        dphi    - interval of azimuth
        group   - compute group velocities or nor
        outfname- output file name (ASDF format)
        --------------------------------------------------------------------------------------------
        Output:
        self.phvelArr       - phace velocity array (3*Nphi*Ntheta)
        self.eigvecArr      - eigenvector array (polarization of three phase modes, 3*3*Nphi*Ntheta)
        self.grvelArr       - group velocity array (3*Nphi*Ntheta)
        self.group_vecArr   - polarization of three group modes (3*3*Nphi*Ntheta)
        self.group_pvArr    - propagation vector of three group modes (3*3*Nphi*Ntheta)
        self.group_thetaArr - group velocity polar angle array (3*Nphi*Ntheta)
        self.group_phiArr   - group velocity azimuth array (3*Nphi*Ntheta)
        self.pf_angleArr    - power flow angle array (3*Nphi*Ntheta)
        ============================================================================================
        """
        Ntheta  = int(180./dtheta) + 1
        Nphi    = int(360./dphi) 
        thetas  = np.linspace(0., 180., num=Ntheta)
        phis    = np.linspace(0., 360., num=Nphi)
        thetaArr, phiArr = np.meshgrid(thetas, phis)
        # initialization
        phvelArr = np.zeros([3, Nphi, Ntheta])
        eigvecArr= np.zeros([3, 3, Nphi, Ntheta])
        if group:
            grvelArr        = np.zeros([3, Nphi, Ntheta])
            group_vecArr    = np.zeros([3, 3, Nphi, Ntheta])
            group_pvArr     = np.zeros([3, 3, Nphi, Ntheta])
            group_thetaArr  = np.zeros([3, Nphi, Ntheta])
            group_phiArr    = np.zeros([3, Nphi, Ntheta])
            pf_angleArr     = np.zeros([3, Nphi, Ntheta])
        print 'Solving Kelvin-Christoffel equation for a sphere of propagation vector'
        N = Nphi*Ntheta
        dN= int(N/10)
        i = 0; j = 0
        for itheta in xrange(Ntheta):
            for iphi in xrange(Nphi):
                i+=1
                if i%dN ==0:
                    j+=10
                    print 'Finished %'+str(j)
                theta       = thetaArr[iphi, itheta]
                phi         = phiArr[iphi, itheta]
                temp_kceq   = self.copy()
                temp_kceq.set_direction_spherical(theta=theta, phi=phi)
                temp_kceq.get_phvel()
                phvelArr[:, iphi, itheta]       = temp_kceq.phvel
                eigvecArr[:, :, iphi, itheta]   = temp_kceq.eig_vec
                if not group: continue
                temp_kceq.get_grvel()
                grvelArr[:, iphi, itheta]       = temp_kceq.grvel
                group_vecArr[:, :, iphi, itheta]= temp_kceq.group_vec
                group_pvArr[:, :, iphi, itheta] = temp_kceq.group_pv
                group_thetaArr[:, iphi, itheta] = temp_kceq.group_theta
                group_phiArr[:, iphi, itheta]   = temp_kceq.group_phi
                pf_angleArr[:, iphi, itheta]    = temp_kceq.powflow_angle
        print 'End solving Kelvin-Christoffel equation for a sphere of propagation vector'
        self.phvelArr   = phvelArr
        self.eigvecArr  = eigvecArr
        if group:
            self.grvelArr        = grvelArr
            self.group_vecArr    = group_vecArr
            self.group_pvArr     = group_pvArr
            self.group_thetaArr  = group_thetaArr
            self.group_phiArr    = group_phiArr
            self.pf_angleArr     = pf_angleArr
        if outfname is not None:
            dset    = pyasdf.ASDFDataSet(outfname)
            dset.add_auxiliary_data(data=thetaArr,  data_type='PySAT', path='theta', parameters={})
            dset.add_auxiliary_data(data=phiArr,  data_type='PySAT', path='phi', parameters={})
            dset.add_auxiliary_data(data=phvelArr,  data_type='PySAT', path='phvel', parameters={})
            dset.add_auxiliary_data(data=eigvecArr,  data_type='PySAT', path='eigvec', parameters={})
            if group:
                dset.add_auxiliary_data(data=grvelArr,  data_type='PySAT', path='grvel', parameters={})
                dset.add_auxiliary_data(data=group_vecArr,  data_type='PySAT', path='grvec', parameters={})
                dset.add_auxiliary_data(data=group_pvArr,  data_type='PySAT', path='grpv', parameters={})
                dset.add_auxiliary_data(data=group_thetaArr,  data_type='PySAT', path='grtheta', parameters={})
                dset.add_auxiliary_data(data=group_phiArr,  data_type='PySAT', path='grphi', parameters={})
                dset.add_auxiliary_data(data=pf_angleArr,  data_type='PySAT', path='pfangle', parameters={})
        return
    
    def read_asdf(self, infname):
        """
        Read sphere data from ASDF file
        """
        dset    = pyasdf.ASDFDataSet(infname)
        self.thetaArr   = dset.auxiliary_data['PySAT']['theta'].data.value
        self.phiArr     = dset.auxiliary_data['PySAT']['phi'].data.value
        self.phvelArr   = dset.auxiliary_data['PySAT']['phvel'].data.value
        self.eigvecArr  = dset.auxiliary_data['PySAT']['eigvec'].data.value
        try:
            self.grvelArr       = dset.auxiliary_data['PySAT']['grvel'].data.value
            self.group_vecArr   = dset.auxiliary_data['PySAT']['grvec'].data.value
            self.group_pvArr    = dset.auxiliary_data['PySAT']['grpv'].data.value
            self.group_thetaArr = dset.auxiliary_data['PySAT']['grtheta'].data.value
            self.group_phiArr   = dset.auxiliary_data['PySAT']['grphi'].data.value
            self.pf_angleArr    = dset.auxiliary_data['PySAT']['pfangle'].data.value
            print 'Group velocity related data available !'
        except AttributeError:
            print 'No group velocity related data !'
        return
    
    def plot3d(self, size=(1000, 700), datatype='phase', ptype='absolute', stype='absolute', fastpolar=True, slowpolar=False, ds=10):
        """
        Plot 3D pole figure using Mayavi 
        ============================================================================================
        Input Parameters:
        size        - size of the Mayavi window
        datatype    - phase or group
        ptype       -   1. 'absolute' or 'abs': plotting P wave in km/s for absolute velocity
                        2. 'relative' or 'rel': plotting P wave in % for anisotropy percentage
        stype       -   1. 'absolute' or 'abs': plotting S wave in km/s for absolute velocity
                        2. 'relative' or 'rel': plotting S wave in % for anisotropy percentage
        fastpolar   - plot polarization of qS1 (fast) wave or not
        slowpolar   - plot polarization of qS2 (slow) wave or not
        ds          - downsampling spacing for polarization vector
        ============================================================================================
        """
        # get data for plot
        # position data
        theta       = self.thetaArr/180.*np.pi
        phi         = self.phiArr/180.*np.pi
        sin_theta   = np.sin(theta)
        cos_theta   = np.cos(theta)
        sin_phi     = np.sin(phi)
        cos_phi     = np.cos(phi)
        x           = cos_phi * sin_theta
        y           = sin_phi * sin_theta
        z           = cos_theta
        xp          = 1.05*cos_phi * sin_theta
        yp          = 1.05*sin_phi * sin_theta
        zp          = 1.05*cos_theta
        if datatype == 'phase':
            s2      = self.phvelArr[0,:,:]
            s1      = self.phvelArr[1,:,:]
            p       = self.phvelArr[2,:,:]
            if fastpolar or slowpolar:
                u2  = self.eigvecArr[0,0,:,:]
                v2  = self.eigvecArr[0,1,:,:]
                w2  = self.eigvecArr[0,2,:,:]
                u1  = self.eigvecArr[1,0,:,:]
                v1  = self.eigvecArr[1,1,:,:]
                w1  = self.eigvecArr[1,2,:,:]
        elif datatype == 'group':
            s2      = self.grvelArr[0,:,:]
            s1      = self.grvelArr[1,:,:]
            p       = self.grvelArr[2,:,:]
            if fastpolar or slowpolar:
                u2  = self.group_vecArr[0,0,:,:]
                v2  = self.group_vecArr[0,1,:,:]
                w2  = self.group_vecArr[0,2,:,:]
                u1  = self.group_vecArr[1,0,:,:]
                v1  = self.group_vecArr[1,1,:,:]
                w1  = self.group_vecArr[1,2,:,:]
        if ds > 1:
            xp      = xp[0:-1:ds, 0:-1:ds]
            yp      = yp[0:-1:ds, 0:-1:ds]
            zp      = zp[0:-1:ds, 0:-1:ds]
            u1      = u1[0:-1:ds, 0:-1:ds]
            v1      = v1[0:-1:ds, 0:-1:ds]
            w1      = w1[0:-1:ds, 0:-1:ds]
            u2      = u2[0:-1:ds, 0:-1:ds]
            v2      = v2[0:-1:ds, 0:-1:ds]
            w2      = w2[0:-1:ds, 0:-1:ds]
        diffs = s1 - s2
        if ptype == 'relative' or ptype == 'rel':
            p   = (p - self.iso_P)/self.iso_P * 100.
        if stype == 'relative'or stype == 'rel':
            s1  = (s1 - self.iso_S)/self.iso_S* 100.
            s2  = (s2 - self.iso_S)/self.iso_S* 100.
            diffs = diffs/self.iso_S* 100.
        #############################
        # qp wave pole figure
        #############################
        mayavi.mlab.figure(figure=None, bgcolor=None, fgcolor=None, engine=None, size=size)
        fig3d=mayavi.mlab.mesh(x, y, z, scalars=p)
        fig3d.module_manager.scalar_lut_manager.reverse_lut = True
        if ptype == 'absolute' or ptype == 'abs':
            cb=mayavi.mlab.colorbar(title=datatype+' velocity (km/s)', orientation='horizontal')
        else:
            cb=mayavi.mlab.colorbar(title=datatype+' velocity anisotropy(%)', orientation='horizontal')
        cb.scalar_bar_representation.proportional_resize=True
        tl=mayavi.mlab.title('qP wave')
        # tl.x_position=0.47
        tl.property.font_size=10
        oaxes=mayavi.mlab.orientation_axes()
        #############################
        # qs1 wave (slow) pole figure
        #############################
        mayavi.mlab.figure(figure=None, bgcolor=None, fgcolor=None, engine=None, size=size)
        fig3d=mayavi.mlab.mesh(x, y, z, scalars=s1)
        fig3d.module_manager.scalar_lut_manager.reverse_lut = True
        if fastpolar:
            mayavi.mlab.quiver3d(xp, yp, zp, u1, v1, w1, line_width=0.01, color=(0, 0, 0), mode='2ddash', scale_factor=0.05)
            mayavi.mlab.quiver3d(xp, yp, zp, -u1, -v1, -w1, line_width=0.01, color=(0, 0, 0), mode='2ddash', scale_factor=0.05)
        if stype == 'absolute' or stype == 'abs':
            cb=mayavi.mlab.colorbar(title=datatype+' velocity (km/s)', orientation='horizontal')
        else:
            cb=mayavi.mlab.colorbar(title=datatype+' velocity anisotropy(%)', orientation='horizontal')
        cb.scalar_bar_representation.proportional_resize=True
        tl=mayavi.mlab.title('qS1(fast) wave',)
        # tl.x_position=0.47
        tl.property.font_size=10
        # mayavi.mlab.text3d(0, 0, 1, 'label')
        oaxes=mayavi.mlab.orientation_axes()
        #############################
        # qs2 wave (fast) pole figure
        #############################
        mayavi.mlab.figure(figure=None, bgcolor=None, fgcolor=None, engine=None, size=size)
        fig3d=mayavi.mlab.mesh(x, y, z, scalars=s2)
        fig3d.module_manager.scalar_lut_manager.reverse_lut = True
        if slowpolar:
            mayavi.mlab.quiver3d(xp, yp, zp, u2, v2, w2, line_width=0.01, color=(1, 1, 1), mode='2ddash', scale_factor=0.05)
            mayavi.mlab.quiver3d(xp, yp, zp, -u2, -v2, -w2, line_width=0.01, color=(1, 1, 1), mode='2ddash', scale_factor=0.05)
        if stype == 'absolute' or stype == 'abs':
            cb=mayavi.mlab.colorbar(title=datatype+' velocity (km/s)', orientation='horizontal')
        else:
            cb=mayavi.mlab.colorbar(title=datatype+' velocity anisotropy(%)', orientation='horizontal')
        cb.scalar_bar_representation.proportional_resize=True
        tl=mayavi.mlab.title('qS2(slow) wave',)
        # tl.x_position=0.47
        tl.property.font_size=10
        oaxes=mayavi.mlab.orientation_axes()
        #############################
        # s wave difference wave pole figure
        #############################
        mayavi.mlab.figure(figure=None, bgcolor=None, fgcolor=None, engine=None, size=size)
        fig3d=mayavi.mlab.mesh(x, y, z, scalars=diffs)
        fig3d.module_manager.scalar_lut_manager.reverse_lut = True
        if fastpolar:
            mayavi.mlab.quiver3d(xp, yp, zp, u1, v1, w1, line_width=0.01, color=(0, 0, 0), mode='2ddash', scale_factor=0.05)
            mayavi.mlab.quiver3d(xp, yp, zp, -u1, -v1, -w1, line_width=0.01, color=(0, 0, 0), mode='2ddash', scale_factor=0.05)
        if slowpolar:
            mayavi.mlab.quiver3d(xp, yp, zp, u2, v2, w2, line_width=0.01, color=(1, 1, 1), mode='2ddash', scale_factor=0.05)
            mayavi.mlab.quiver3d(xp, yp, zp, -u2, -v2, -w2, line_width=0.01, color=(1, 1, 1), mode='2ddash', scale_factor=0.05)
        if stype == 'absolute' or stype == 'abs':
            cb=mayavi.mlab.colorbar(title=datatype+' velocity (km/s)', orientation='horizontal')
        else:
            cb=mayavi.mlab.colorbar(title=datatype+' velocity anisotropy(%)', orientation='horizontal')
        cb.scalar_bar_representation.proportional_resize=True
        tl=mayavi.mlab.title('S wave difference',)
        # tl.x_position=0.47
        tl.property.font_size=10
        oaxes=mayavi.mlab.orientation_axes()
        return
    
    def plot2d(self, size=(10,10), cmap='jet_r', hsph='upper', contour=False, theta0=180., phi0=0., datatype='phase', ptype='absolute', stype='absolute', fastpolar=True, slowpolar=False, ds=10):
        """
        Plot 2D pole figure using cartopy 
        ============================================================================================
        Input Parameters:
        size        - size of the figure
        cmap        - colormap type
        contour     - plot contour or not
        theta0, phi0- center polar angle and azimuth for view
        datatype    - phase or group
        ptype       -   1. 'absolute' or 'abs': plotting P wave in km/s for absolute velocity
                        2. 'relative' or 'rel': plotting P wave in % for anisotropy percentage
        stype       -   1. 'absolute' or 'abs': plotting S wave in km/s for absolute velocity
                        2. 'relative' or 'rel': plotting S wave in % for anisotropy percentage
        fastpolar   - plot polarization of qS1 (fast) wave or not
        slowpolar   - plot polarization of qS2 (slow) wave or not
        ds          - downsampling spacing for polarization vector
        ============================================================================================
        """
        # Get data for plot
        # Position data
        import pycpt
        if cmap == 'cv':
            try: cmap=pycpt.load.gmtColormap('cv.cpt')
            except: cmap='jet_r'
        if cmap == 'lasif':
            from lasif import colors
            cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        lonArr = self.phiArr.T; latArr = (90.-self.thetaArr).T
        phip   = self.phiArr.T; thetap = self.thetaArr.T
        lonp = lonArr.copy(); latp=latArr.copy()
        lon0=phi0; lat0=90.-theta0
        if datatype == 'phase':
            s2      = (self.phvelArr[0,:,:]).T
            s1      = (self.phvelArr[1,:,:]).T
            p       = (self.phvelArr[2,:,:]).T
            if fastpolar or slowpolar:
                u2  = (self.eigvecArr[0,0,:,:]).T
                v2  = (self.eigvecArr[0,1,:,:]).T
                w2  = (self.eigvecArr[0,2,:,:]).T
                u1  = (self.eigvecArr[1,0,:,:]).T
                v1  = (self.eigvecArr[1,1,:,:]).T
                w1  = (self.eigvecArr[1,2,:,:]).T
        elif datatype == 'group':
            s2      = (self.grvelArr[0,:,:]).T
            s1      = (self.grvelArr[1,:,:]).T
            p       = (self.grvelArr[2,:,:]).T
            if fastpolar or slowpolar:
                u2  = (self.group_vecArr[0,0,:,:]).T
                v2  = (self.group_vecArr[0,1,:,:]).T
                w2  = (self.group_vecArr[0,2,:,:]).T
                u1  = (self.group_vecArr[1,0,:,:]).T
                v1  = (self.group_vecArr[1,1,:,:]).T
                w1  = (self.group_vecArr[1,2,:,:]).T
        if ds > 1:
            lonp    = lonp[0:-1:ds, 0:-1:ds]
            latp    = latp[0:-1:ds, 0:-1:ds]
            phip    = phip[0:-1:ds, 0:-1:ds]
            thetap  = thetap[0:-1:ds, 0:-1:ds]
            if fastpolar or slowpolar:
                u1      = u1[0:-1:ds, 0:-1:ds]
                v1      = v1[0:-1:ds, 0:-1:ds]
                w1      = w1[0:-1:ds, 0:-1:ds]
                u2      = u2[0:-1:ds, 0:-1:ds]
                v2      = v2[0:-1:ds, 0:-1:ds]
                w2      = w2[0:-1:ds, 0:-1:ds]
        diffs = s1 - s2
        if ptype == 'relative' or ptype == 'rel':
            p   = (p - self.iso_P)/self.iso_P * 100.
        if stype == 'relative'or stype == 'rel':
            s1  = (s1 - self.iso_S)/self.iso_S* 100.
            s2  = (s2 - self.iso_S)/self.iso_S* 100.
            diffs = diffs/self.iso_S* 100.
        if fastpolar:
            rr1, uu1, vv1 =  localcartesian2spherical(u1, v1, w1, thetap, phip, r=100)
        if slowpolar:
            rr2, uu2, vv2 =  localcartesian2spherical(u2, v2, w2, thetap, phip, r=100)
        if hsph == 'lower':
            p   = p[::-1, :]
            s1  = s1[::-1, :]
            s2  = s2[::-1,:]
            print fastpolar,slowpolar
        
            if fastpolar:
              uu1 = uu1[::-1, :]
              vv1 = vv1[::-1, :]
            if slowpolar:
              uu2 = uu2[::-1,:]
              vv2 = vv2[::-1, :]
              
        fig = plt.figure(figsize=size)
        #############################
        # qP wave pole figure
        #############################
        plt.subplot(221, projection=ccrs.Orthographic(central_longitude=lon0, central_latitude=lat0))
        if not contour:
            im=plt.pcolormesh(lonArr, latArr, p, transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            levels=np.linspace(p.min(), p.max(), 20)
            im=plt.contourf(lonArr, latArr, p,  transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
        # return lonp, latp, uu1, vv1, rr1
        cb = plt.colorbar(im)
        if ptype == 'absolute' or ptype == 'abs':
            cb.set_label(datatype+' velocity (km/s)', fontsize=15, rotation=90)
        else:
            cb.set_label(datatype+' velocity anisotropy(%)', fontsize=15, rotation=90)
        plt.title('qP wave', fontsize=20)
        #############################
        # qS1 wave pole figure
        #############################
        plt.subplot(222, projection=ccrs.Orthographic(central_longitude=lon0, central_latitude=lat0))
        if not contour:
            im=plt.pcolormesh(lonArr, latArr, s1, transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            levels=np.linspace(p.min(), p.max(), 20)
            im=plt.contourf(lonArr, latArr, s1,  transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
        cb = plt.colorbar(im)
        if stype == 'absolute' or stype == 'abs':
            cb.set_label(datatype+' velocity (km/s)', fontsize=15, rotation=90)
        else:
            cb.set_label(datatype+' velocity anisotropy(%)', fontsize=15, rotation=90)
        if fastpolar:
            plt.quiver(lonp, latp, uu1, vv1, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(0,0,0))
            plt.quiver(lonp, latp, -uu1, -vv1, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(0,0,0))
        plt.title('qS1(fast) wave', fontsize=20)
        #############################
        # qS2 wave pole figure
        #############################
        plt.subplot(223, projection=ccrs.Orthographic(central_longitude=lon0, central_latitude=lat0))
        if not contour:
            im=plt.pcolormesh(lonArr, latArr, s2, transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            levels=np.linspace(p.min(), p.max(), 20)
            im=plt.contourf(lonArr, latArr, s2,  transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
        cb = plt.colorbar(im)
        if stype == 'absolute' or stype == 'abs':
            cb.set_label(datatype+' velocity (km/s)', fontsize=15, rotation=90)
        else:
            cb.set_label(datatype+' velocity anisotropy(%)', fontsize=15, rotation=90)
        if slowpolar:
            plt.quiver(lonp, latp, uu2, vv2, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(1,1,1))
            plt.quiver(lonp, latp, -uu2, -vv2, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(1,1,1))
        plt.title('qS2(slow) wave', fontsize=20)
        #############################
        # S wave difference pole figure
        #############################
        plt.subplot(224, projection=ccrs.Orthographic(central_longitude=lon0, central_latitude=lat0))
        if not contour:
            im=plt.pcolormesh(lonArr, latArr, diffs, transform=ccrs.PlateCarree(), cmap=cmap)
        else:
            levels=np.linspace(p.min(), p.max(), 20)
            im=plt.contourf(lonArr, latArr, diffs,  transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
        cb = plt.colorbar(im)
        if stype == 'absolute' or stype == 'abs':
            cb.set_label(datatype+' velocity (km/s)', fontsize=15, rotation=90)
        else:
            cb.set_label(datatype+' velocity anisotropy(%)', fontsize=15, rotation=90)
        if fastpolar:
            plt.quiver(lonp, latp, uu1, vv1, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(0,0,0))
            plt.quiver(lonp, latp, -uu1, -vv1, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(0,0,0))
        if slowpolar:
            plt.quiver(lonp, latp, uu2, vv2, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(1,1,1))
            plt.quiver(lonp, latp, -uu2, -vv2, transform=ccrs.PlateCarree(), scale=50, width=0.01, headaxislength=0, headlength=0, color=(1,1,1))
        plt.title('S wave difference', fontsize=20)
        plt.suptitle(self.etensor.info, fontsize=15)
        plt.show()
        return
        


