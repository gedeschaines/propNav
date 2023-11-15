# -*- coding: utf-8 -*-

# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name
# pylint: disable=anomalous-backslash-in-string,bad-continuation
# pylint: disable=multiple-statements,redefined-outer-name,global-statement

# File: propNav.py
# Auth: Gary E. Deschaines
# Date: 25 Oct 2023
# Prog: Proportional navigation guidance missile flyout model
# Desc: Application of selectable proportional navigation guidance laws
#       (True, Pure, ZEM, or Augmented PN) for missile engagement
#       of a target. 3-DOF point mass kinematic model for missile and
#       target. Ideal missile control; no lag with 100% effective, but 
#       bounded commanded acceleration, and perfect command response.
#
# Note: Refactored from a Mathcad 3-DOF kinematic ideal proportional
#       navigation guidance missile flyout model developed in 1997.
#       The inertial (fixed) reference frame Cartesian (+X, +Y, +Z)
#       coordinate system in the Mathcad model correlates with (East,
#       North, Up), while translational/rotational missile and target
#       body frame (+x, +y, +z) axes follow the (forward, right, down)
#       convention. Care must be taken in transforming displacement 
#       and directional vectors between these coordinate frames, and
#       describing rotational directions. Specifically, since missile
#       and target body frame +x axes are aligned with their respective
#       inertial velocity vectors, positive azimuth rotation is negative
#       body yaw, while positive elevation rotation is positive body
#       pitch. Thus, positive accelerations normal to body frame +x
#       axis are those resulting from positive yaw or pitch rates 
#       crossed with the body frame inertial velocity vector. The +/-
#       missile line-of-sigt (LOS) rates and normal accelerations will
#       be evident in associated profile plots.
#
# Refs:
#
#  [1] Paul Zarchan and A. Richard Seebass (Editor-in-Chief), 
#      "Tactical and Strategic Missile Guidance (Progress in 
#      Astronautics and Aeronautics, Vol 124)", American 
#      Institute of Aeronautics and Astronautics, Washington,
#      D.C., 1990.
#
#  [2] Donald T. Greenwood, "Principles of Dynamics", Prentice-Hall,
#      Inc. of Englewood Clifts, New Jersey, 1965.
#
#  [3] Neil F. Palumbo, Ross A. Blauwkamp, and Justin M. Lloyd,
#      "Basic Principles of Homing Guidance", rev 2018, Johns
#      Hopkins APL Technical Digest, VOL 29, No 1, 2010. Web
#      available at secwww.jhuapl.edu/techdigest:
#      https://secwww.jhuapl.edu/techdigest/Content/techdigest/pdf/V29-N01/29-01-Palumbo_Principles_Rev2018.pdf)
#
#  [4] Ben Dickinson, "Missile Guidance Fundamentals Tutorial", 
#      last updated Oct. 15, 2023. Web available at www.youtube.com:
#      https://www.youtube.com/playlist?list=PLcmbTy9X3gXt02z1wNy4KF5ui0tKxdQm7)
#
#  [5] Ben Dickinson, "Guidance from Optimal Control",
#      last updated Apr. 2, 2023. Web available at www.youtube.com:
#      https://www.youtube.com/playlist?list=PLcmbTy9X3gXsh-o1W60E7rEA35igNj__q)
#
#  [6] Farham A. Faruqi, "Integrated Navigation, Guidance, and
#      Control of Missile Systems: 3-D Dynamic Model", Weapon 
#      Systems Division DSTO, DSTO-TR-2805, Feb., 2013. Web 
#      available at www.dst.defence.gov.au:
#      https://www.dst.defence.gov.au/publication/integrated-navigation-guidance-and-control-missile-systems-3-d-dynamic-model
#
# Disclaimer:
#
# See DISCLAIMER

import sys


from math import ceil, floor, cos, sin, acos, asin, atan, atan2, pi 
from io import StringIO
#from locale import format_string

try:
    import numpy             as np
    import numpy.linalg      as la
    import matplotlib        as mpl
    import matplotlib.pyplot as plt
#   import matplotlib.animation as animation
#   from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d.art3d import Line3D
except ImportError:
    print("* Error: NumPy and Matplotlib packages required.")
    sys.exit()

try:
    from RK4_Solver import RK4_Solver
except ImportError:
    print("* Error: RK4_Solver class required.")
    sys.exit()


###
### Global constants
###

RPD = atan(1.0)/45.0  # radians per degree
DPR = 1.0/RPD         # degrees per radian
g   = 9.80665         # gravitational acceleration at sea-level (meters/s/s)

# Unit vectors for inertial Cartesion frame X,Y,Z axes.

Uxi = np.array([1.0, 0.0, 0.0])
Uyi = np.array([0.0, 1.0, 0.0])
Uzi = np.array([0.0, 0.0, 1.0])  

# Proportional Navigation law (method) selection.

PN_TRUE = 1  # Default
PN_PURE = 2
PN_ZEM  = 3  # Zero Effort Miss
PN_ATPN = 4  # Augmented True Proportional Navigation
PN_APPN = 5  # Augmented Pure Proportional Navigation
PN_LAWS = {PN_TRUE:'True', PN_PURE:'Pure', PN_ZEM:'ZEM', 
           PN_ATPN:'ATPN', PN_APPN:'APPN'}
PNAV    = PN_PURE

Nm = 4    # proportional navigation constant
Nt = 3.0  # target turning acceleration (g's)

# Set missile type and acceleration maximum.

SAM = 1  # For engagements described in Caveats section of propNav README.
AAM = 2  # For engagements presented in Section 3, Modules 3 & 4, Section 4,
         # Module 4 of ref [4], and Section 2, Module 3 of ref [5].
MSL = SAM

Gmmax = {SAM:8, AAM:30}  # maximum missile acceleration (g's)
Ammax = Gmmax[MSL]*g     # maximum missile acceleration (meters/s/s)

# Set minimum miss distance (meters).

if MSL == SAM:
    MinMissDist = 3.0
else :
    MinMissDist = 6.0
    
# Missile lead azimuth and elevation angles in degrees. 
# Set to None for calculation of lead angles based on
# estimation of time-to-intercept of a non-maneuvering
# target with constant velocity and heading.

if MSL == SAM:
   maz = 10.0
   mel = 12.0
else:
   maz = 0.0
   mel = 0.0

# Define target and missile initial states.

if MSL == SAM:
    Pt0   = np.array([ 2000.0,    0.0,  500.0])
    Vt0   = np.array([    0.0,  200.0,    0.0])
    Pm0   = np.array([    0.0,    0.0,    2.0])
    magVm = 450.0
else:
    if ((int(Nt) == 3) and (int(Nm) == 3)) and \
       ((PNAV == PN_TRUE) or (PNAV == PN_ATPN) or (PNAV == PN_APPN)):
        # for Section 2, Module 3 of ref [5].
        Pt0   = np.array([ 9144.0, 4572.0,    0.0])
        Vt0   = np.array([ -304.8,    0.0,    0.0])
        Pm0   = np.array([    0.0, 4572.0,    0.0])
        magVm = 457.2
    else:
        # For Section 3, Modules 3 & 4, Section 4, Module 4 of ref [4].
        Pt0   = np.array([12192.0, 6096.0, 3048.0])
        Vt0   = np.array([ -304.8,    0.0,    0.0])
        Pm0   = np.array([    0.0, 6096.0, 3048.0])
        magVm = 914.4

# Define target turning/climbing rotation axis unit vector.

if MSL == SAM:
    UWt = Uzi  # for level turn in XY plane
    #UWt = np.array([-0.2418, 0.2418, 0.9397])  # for diving turn
else:
    UWt = Uyi  # for climbing turn in XZ plane
UWt = UWt/la.norm(UWt)

# Set integration time step size and simulation stop time (sec).

T_STEP = 0.001
if MSL == SAM:
    T_STOP =  5.5
else:
    if ((int(Nt) == 3) and (int(Nm) == 3)) and \
        ((PNAV == PN_TRUE) or (PNAV == PN_ATPN) or (PNAV == PN_APPN)):
        # For Section 2, Module 3 of ref [5].
        T_STOP = 14.5
    else:
        # For Section 3, Modules 3 & 4, Section 4, Module 4 of ref [4].
        T_STOP = 12.5
    
# Set processing output control flags.

PRINT_DATA = False   # Controls printing of collected data
PLOT_DATA  = True    # Controls plotting of collected data
SAVE_ANIM  = False   # Controls saving animation images (not implemented)
PRINT_TXYZ = True    # Controls printing TXYZ.OUT file

###
### Procedures for Proportional Navigation Guidance model
###

def Uvec(V):
    magV = la.norm(V)
    if magV > 0.0: U = V / magV
    else:          U = np.array([0.0, 0.0, 0.0])
    return U

def Prel(Pt, Pm):
    return Pt - Pm

def Vrel(Vt, Vm):
    return Vt - Vm

def Mrot(psi, tht, phi):
    """
    Euler angle (yaw,pitch,roll) rotation transformation matrix
    for inertial frame Cartesian coordinates (Xi,Yi,Zi) to body
    frame coordinates (xb,yb,zb) derived from equation (7-166)
    on pg 335 of ref [2] (same as that presented in Figure 3 on
    pg 6 of ref [6]). In the derivation below, consider the 
    xyz inertial frame depicted in Figure 3 on pg 6 of ref [6]
    rotated 180 deg counter-clockwise about the +yi axis so the
    +zi axis points toward top of the page, then the following
    relationships exist between azimuth, elevation and bank 
    angles in propNav inertial frame and psi, theta, phi in 
    the figure:
        
      psi = azimuth
      theta = -elevation
      phi = bank angle
     
    Note:  Yaw, pitch and roll are expected to be in radians,
           NOT degrees!
    
    Use as:  [xb, yb, zb] = Numpy.matmul(Mbi, [Xi, Yi, Zi])
    """
    cpsi = cos(psi)
    spsi = sin(psi)
    ctht = cos(-tht)
    stht = sin(-tht)
    cphi = cos(phi)
    sphi = sin(phi)
    
    M = np.zeros([3,3])
    
    M[0,0] =  cpsi*ctht
    M[0,1] =  spsi*ctht
    M[0,2] = -stht
    M[1,0] = -spsi*cphi + cpsi*stht*sphi
    M[1,1] =  cpsi*cphi + spsi*stht*sphi
    M[1,2] =  ctht*sphi
    M[2,0] =  spsi*sphi + cpsi*stht*cphi
    M[2,1] = -cpsi*sphi + spsi*stht*cphi
    M[2,2] =  ctht*cphi  
    
    # Account for +yb = -Yi and +zb = -Zi
    Mbi = np.matmul(np.array([[1.0,  0.0,  0.0],
                              [0.0, -1.0,  0.0],
                              [0.0,  0.0, -1.0]]), M)
    return Mbi

def leadAngle(Pt, Vt, Pm, magVm):
    #
    # Assuming constant velocity magnitudes and
    # headings, at some future intercept point the
    # triangle formed by missile velocity magVm and
    # target velocity magVt toward that point from 
    # their current positions Pm and Pt, and the
    # LOS distance between current positions (i.e., 
    # Rlos = Pt - Pm) must satisfy the following
    # relationship:
    #
    #   Vm*sin(alpha) = Vt*sin(beta)
    #
    # where alpha (lead angle) is the angle between 
    # Vm and Rlos, and beta is the angle between Vt
    # and Rlos. Using definition of inner product:
    #
    #   cos(beta) = <Vt, Rlos> / |Vt||Rlos|
    #
    # a value for beta can determined, and then the
    # above relationship solved for alpha.
    #
    magVt = la.norm(Vt)
    Ptm   = Prel(Pt, Pm)
    cosb  = np.dot(Vt, Ptm)/(magVt*la.norm(Ptm))
    alpha = asin((magVt/magVm)*sin(acos(cosb)))
    return alpha  # Note: alpha in radians.

def az_el_of_V(V):
    # Note: DPR is global
    U  = Uvec(V)
    az = atan2(U[1], U[0])*DPR
    el = atan(U[2]/la.norm([U[0], U[1]]))*DPR
    return az, el  # Note: az and el in degrees.

def setVm(vmag, az, el):
    vx = vmag*cos(el)*cos(az)
    vy = vmag*cos(el)*sin(az)
    vz = vmag*sin(el)
    Vm = np.array([vx, vy, vz])
    return Vm

def Vclose(Vt, Vm, Ulos, collision=False):
    # Note:  Closing velocity is defined as -d(Rlos)/dt.
    if collision == True:
        # Calculate collision closing velocity (see calcVcTgo).
        vt = la.norm(Vt)
        vm = la.norm(Vm)
        betat = acos(np.dot(Vt, Ulos)/vt)
        betam = asin((vt/vm)*sin(betat))
        Vc = (vm*cos(betam) - vt*cos(betat))*Ulos
    else:
        # Calculate standard closing velocity which assumes 
        # Vt, Vm and Ulos are all within the same plane.
        Vc = -np.dot(Vrel(Vt,Vm), Ulos)*Ulos
    return Vc
 
def timeToGo(Rlos, Vt, Vm):
    # Note:  Uses collision closing velocity (see calcVcTgo).
    vcc = la.norm(Vclose(Vt, Vm, Uvec(Rlos), True))
    tgo = la.norm(Rlos)/vcc
    return tgo

def calcVcTgo(Pt, Vt, Pm, Vm):
    #
    # Calculates Vc as collision course relative (closing)
    # velocity of missile wrt target, and Tgo as time-to-go
    # (to intercept) using equations presented on pgs 25-26
    # in sections C.1.1, C.1.2 and C.1.3 of ref [6].
    #
    azt, elt   = az_el_of_V(Vt)
    aztm, eltm = az_el_of_V(Prel(Pt,Pm))
    
    thtt  = elt*RPD
    psit  = -azt*RPD
    thttm = eltm*RPD
    psitm = -aztm*RPD
    
    # Unit vector (evt) along target body, and unit Rlos (estm) vector;
    # equations (A2.1) and (A2.2).
    evt  = np.array([cos(thtt)*cos(psit), cos(thtt)*sin(psit), -sin(thtt)])
    estm = np.array([cos(thttm)*cos(psitm), cos(thttm)*sin(psitm), -sin(thttm)])
 
    # Angle between target velocity vector Vt and the Rlos vector measured
    # in Vt X Rlos X Vm plane; equations (A2.3) and (A2.4).
    Betatm = acos(np.dot(evt, estm))
    
    # Angle between missile collision course velocity vector and the
    # RLos vector measured in Vt X Rlos X Vm plane; equation (A2.6).
    Betaccmt = asin((la.norm(Vt)/la.norm(Vm))*sin(Betatm))
    
    # Collision course closing velocity of missile wrt target along
    # Rmt; equation (A2.7).
    VCccmt = la.norm(Vm)*cos(Betaccmt) - la.norm(Vt)*cos(Betatm)
    
    # Target/missile range-to-go; equation (A2.8).
    Rmt = la.norm(Prel(Pm, Pt))
    
    # Time-to-go; equation (A2.9).
    Tgo = Rmt/VCccmt
    
    return VCccmt, Tgo

def ZEMn(Rlos, Vtm, tgo):
    Ulos = Uvec(Rlos)
    ZEM  = Rlos + Vtm*tgo
    ZEMr = np.dot(ZEM, Ulos)*Ulos
    ZEMn = ZEM - ZEMr 
    return ZEMn

def Wlos(Vt, Vm, Rlos, Ulos):
    # Note:  Following expressions for calculating Wlos
    #
    Vnrm = Vrel(Vt, Vm) - Vclose(Vt, Vm, Ulos)
    Wlos = np.cross(-Vnrm, Ulos)/la.norm(Rlos)
    #
    # is equivalent to:
    #
    #   Wlos = np.cross(Rlos, Vtm)/np.dot(Rlos, Rlos)
    #
    # which can be reduced to:
    #
    #   Wlos = np.array([Rlos[1]*Vtm[2] - Rlos[2]*Vtm[1],
    #                    Rlos[2]*Vtm[0] - Rlos[0]*Vtm[2],
    #                    Rlos[0]*Vtm[1] - Rlos[1]*Vtm[0]
    #                   ])/np.dot(Rlos, Rlos)
    #
    # as shown in derivation of equation (2.18) in ref [6].
    return Wlos

def Amslc(Rlos, Vt, At, Vm, N):
    """
    This routine is the application of selected proportional
    navigation method - True, Pure, FEM, ATPN or APPN.
    
    Globals
    -------
    PNAV : integer constant 
        Proportional Navigation law selected identifier
    RPD : float constant
        Radians per degree
        
    Parameters
    ----------
    Rlos : float 3-vector
        Range along LOS from missile to target.
    Vt : float 3-vector
        Velocity (inertial) of target.
    At : float 3-vector
        Acceleration (inertial) of target.
    Vm : float 3-vector
        Velocity (inertial) of missile.
    N : float
        Proportional navigation constant (or gain).

    Returns
    -------
    Acmd : float 3-Vector
        Acceleration commanded.
        
    """
    Ulos = Uvec(Rlos)
    Vtm  = Vrel(Vt, Vm)
    if PNAV == PN_APPN or PNAV == PN_ATPN:
        # Create inertial to missile body rotation matrix
        maz, mel = az_el_of_V(Vm)
        Mbi = Mrot(maz*RPD, mel*RPD, 0.0)
        # See derivation of equation (3.8) in ref [6].
        #
        if PNAV == PN_APPN:
            # 3.1.1 Version 1 (PN-1) Pure PN equation (3.2)
            PN_1i = N*np.cross(Wlos(Vt, Vm, Rlos, Ulos), Vm)
            PN_1b = np.matmul(Mbi, PN_1i)  # eq. (3.3)
            PN_1b[0] = 0.0                 # eq. (3.4)
            PNG = PN_1b
        if PNAV == PN_ATPN:
            # 3.1.2 Version 2 (PN-2) True PN equation (3.5)
            Vc    = la.norm(Vclose(Vt, Vm, Ulos))
            PN_2i = N*Vc*np.cross(Wlos(Vt, Vm, Rlos, Ulos), Ulos)
            PN_2b = np.matmul(Mbi, PN_2i)  # eq. (3.6)
            PN_2b[0] = 0.0                 # eq. (3.7)
            PNG = PN_2b  # Section 2, Module 3 of ref [5].
        #
        # 3.2 Augmented PN (APN) Guidance
        Acmdb = PNG + np.matmul(Mbi, (N/2.0)*At)  # eq. (3.8)
        Acmdb[0] = 0.0  # no thrust control
        Acmd = np.matmul(Mbi.transpose(), Acmdb)
    elif PNAV == PN_ZEM:
        tgo  = timeToGo(Rlos, Vt, Vm)
        Acmd = N*ZEMn(Rlos, Vtm, tgo)/(tgo**2)
        Acmd = Acmd - np.dot(Acmd, Uvec(Vm))*Uvec(Vm)  # no thrust control
    elif PNAV == PN_PURE:
        Acmd = N*np.cross(Wlos(Vt, Vm, Rlos, Ulos), Vm)
    else: # PNAV == PN_TRUE
        Vc   = la.norm(Vclose(Vt, Vm, Ulos))
        Acmd = N*Vc*np.cross(Wlos(Vt, Vm, Rlos, Ulos), Ulos)
        Acmd = Acmd - np.dot(Acmd, Uvec(Vm))*Uvec(Vm)  # no thrust control
    return Acmd

def Amsla(Amcmd, Ammax):
    """
    Applies Ammax bound to given commanded missile acceleration.

    Parameters
    ----------
    Amcmd : float 3-vector
        Missile inertial linear acceleration commanded.
    Ammax : float constant
        Maximum missile linear acceleration.

    Returns
    -------
    float 3-Vector
        Missile inertial acceleration achieved (actual).

    """
    if la.norm(Amcmd) > Ammax:
        return np.dot(Ammax, Uvec(Amcmd))
    return Amcmd

def Atgt(UWt, Vt, n):
    """
    This routine calculates target inertial linear acceleration
    for given inertial angular velocity direction unit vector,
    inertial linear velocity, and turning g's.
    
    Globals
    -------
    g : float
        gravitional acceleration magnitude

    Parameters
    ----------
    UWt : float 3-vector
        Target angular velocity direction unit vector (i.e.,
        direction frame rotation axis points in inertial space).
    Vt : float 3-vector
        Target inertial velocity.
    n : float 
        Target turning acceleration (normal to UWt) magnitude 
        in g's.

    Returns
    -------
    At : float 3-vector
        Target inertial acceleration.

    """
    if n != 0.0:
        magVt = la.norm(Vt)
        # Calculate lift g's loss due to pitch angle
        theta = pitchAngle(Vt)
        loss  = sin(theta)
        # Calculate target az and el in degrees
        taz, tel = az_el_of_V(Vt)
        # Set target body frame rotation axis
        if abs(np.dot(UWt,Uyi)) < 0.0001:
            UWtb = -UWt  # +yawing about -zb axis
        else:
            UWtb = UWt   # +pitching about +zb axis 
        # Rotate inertial Vt into target body frame
        Mbi = Mrot(taz*RPD, tel*RPD, 0.0)
        Vtb = np.matmul(Mbi, Vt)
        # Calculate inertial acceleration in target body frame
        Atb = np.cross((((n-loss)*g)/magVt)*UWtb, Vtb.flatten())
        # Rotate inertial accelaration into inertial space frame
        At = np.matmul(Mbi.transpose(), Atb)
    else:
        At = np.array([0.0, 0.0, 0.0])
    return At

def pitchAngle(Vt):   
    theta = atan2(Vt[2], la.norm([Vt[0], Vt[1]]))
    return theta  # Note: theta in radians.

def bankAngle(At, Vt):
    # Note: RPD and g are global.
    UVt = Uvec(Vt)
    # Calculate turning acceleration normal to Vt.
    Atn = At - np.dot(At, UVt)*UVt
    # Calculate target az and el in degrees
    taz, tel = az_el_of_V(Vt)
    # Rotate inertial Atn into target body frame
    Mbi = Mrot(taz*RPD, tel*RPD, 0.0)
    Atb = np.matmul(Mbi, Atn)
    # Only use Y component of turning acceleration normal.
    phi = atan(Atb[1]/g)
    return phi  # Note: phi in radians.


###
### Procedures for differential equation integration
###

# Initial values for state variables array.

nSvar = 13               # number of state variables
S     = np.zeros(nSvar)  # state variables
dS    = np.zeros(nSvar)  # state derivatives

def setS(S, Vt, Pt, Vm, Pm):
    # Note: S[0] is t (time)
    S[1]  = Vt[0]
    S[2]  = Vt[1]
    S[3]  = Vt[2]
    S[4]  = Pt[0]
    S[5]  = Pt[1]
    S[6]  = Pt[2]
    S[7]  = Vm[0]
    S[8]  = Vm[1]
    S[9]  = Vm[2]
    S[10] = Pm[0]
    S[11] = Pm[1]
    S[12] = Pm[2]
    return S

def getVtOfS(S):
    return np.array([S[1], S[2], S[3]])
def getPtOfS(S):
    return np.array([S[4], S[5], S[6]])
def getVmOfS(S):
    return np.array([S[7], S[8], S[9]])
def getPmOfS(S):
    return np.array([S[10], S[11], S[12]])

def getS(S):
    Vt = getVtOfS(S)
    Pt = getPtOfS(S)
    Vm = getVmOfS(S)
    Pm = getPmOfS(S)
    return Vt, Pt, Vm, Pm

def setSdot(Sdot, At, Vt, Am, Vm):
    # Note: Sdot[0] is dt (1.0)
    Sdot[1]  = At[0]
    Sdot[2]  = At[1]
    Sdot[3]  = At[2]
    Sdot[4]  = Vt[0]
    Sdot[5]  = Vt[1]
    Sdot[6]  = Vt[2]
    Sdot[7]  = Am[0]
    Sdot[8]  = Am[1]
    Sdot[9]  = Am[2]
    Sdot[10] = Vm[0]
    Sdot[11] = Vm[1]
    Sdot[12] = Vm[2]
    return Sdot

def getAtOfSdot(Sdot):
    return np.array([Sdot[1], Sdot[2], Sdot[3]])
def getVtOfSdot(Sdot):
    return np.array([Sdot[4], Sdot[5], Sdot[6]])
def getAmOfSdot(Sdot):
    return np.array([Sdot[7], Sdot[8], Sdot[9]])
def getVmOfSdot(Sdot):
    return np.array([Sdot[10], Sdot[11], Sdot[12]])

def getSdot(Sdot):
    At = getAtOfSdot(Sdot)
    Vt = getVtOfSdot(Sdot)
    Am = getAmOfSdot(Sdot)
    Vm = getVmOfSdot(Sdot)
    return At, Vt, Am, Vm

def dotS2(dS, At, Vt, Amcmd, Vm):
    # Note: Ammax is global
    dS = setSdot(dS, At, Vt, Amsla(Amcmd, Ammax), Vm)
    return dS

def dotS1(dS, Ptm, Vt, Vm):
    # Note: UWt, Nt and Nm are global
    At = Atgt(UWt ,Vt, Nt)
    dS = dotS2(dS, At, Vt, Amslc(Ptm, Vt, At, Vm, Nm), Vm)
    return dS

def dotS(n, S):
    """
    State derivatives function.
    
    Parameters
    ----------
    n : integer
        Number of state variables (i.e., length of state vector S).
    S : float n-vector
        State variable vector.
        
    Returns
    -------
    dS : float n-Vector
         Computed derivatives of state vector variables (i.e, dS/dt).
         
    """
    dS    = np.zeros(n)
    dS[0] = 1.0  # d(t)/dt
    
    Vt, Pt, Vm, Pm = getS(S)

    dS = dotS1(dS, Prel(Pt,Pm), Vt, Vm)
    
    return dS

def Stop(S):
    Vt, Pt, Vm, Pm = getS(S)
    if (S[0] >= T_STOP) or (np.dot(Prel(Vm,Vt), Uvec(Prel(Pt, Pm))) < 0.0):
        return True
    return False

def Dclose(S):
    Vt, Pt, Vm, Pm = getS(S)
    dtm = la.norm(Prel(Pt, Pm))
    return dtm

def delT(S):
    if Dclose(S) < 100.0:
        if Dclose(S) < 10.0:
            h = 0.00005
        else:
            h = 0.0001
    else:
        h = 0.0005  # Note: this value should be tStep/2.
    return h


# Differential equations of motion constants.
    
tStep = T_STEP  # simulation and integration time step (sec)
tStop = T_STOP  # simulation and integration stop time (sec)

# Animation display and data sampling constants.

FPS    = 100
F_TIME = 1.0/FPS
N_STEP = int(floor((F_TIME + 0.5*T_STEP)/T_STEP))
N_TIME = N_STEP*T_STEP

# Create simulation data collection arrays for plotting.

nSamples = int(ceil(T_STOP/N_TIME)) + 1
    
if PLOT_DATA or PRINT_TXYZ:
    Time  = np.zeros(nSamples+1)  # simulation time
    Ptx   = np.zeros(nSamples+1)  # target position x
    Pty   = np.zeros(nSamples+1)  # target position y
    Ptz   = np.zeros(nSamples+1)  # target position z
    Vtx   = np.zeros(nSamples+1)  # target velocity x
    Vty   = np.zeros(nSamples+1)  # target velocity y
    Vtz   = np.zeros(nSamples+1)  # target velocity z
    Atx   = np.zeros(nSamples+1)  # target acceleration x
    Aty   = np.zeros(nSamples+1)  # target acceleration y
    Atz   = np.zeros(nSamples+1)  # target acceleration z
    Pmx   = np.zeros(nSamples+1)  # missile position x
    Pmy   = np.zeros(nSamples+1)  # missile position y
    Pmz   = np.zeros(nSamples+1)  # missile position z
    Vmx   = np.zeros(nSamples+1)  # missile velocity x
    Vmy   = np.zeros(nSamples+1)  # missile velocity y
    Vmz   = np.zeros(nSamples+1)  # missile velocity z
    Amx   = np.zeros(nSamples+1)  # missile acceleration x
    Amy   = np.zeros(nSamples+1)  # missile acceleration y
    Amz   = np.zeros(nSamples+1)  # missile acceleration z
    Dcls  = np.zeros(nSamples+1)  # Closest approach distance
    LOSd  = np.zeros(nSamples+1)  # LOS rate in (deg/sec)
    VELc  = np.zeros(nSamples+1)  # Closing velocity (meters/sec)
    Acmg  = np.zeros(nSamples+1)  # Missile acceleration in g's
    Velm  = np.zeros(nSamples+1)  # Missile velocity magnitude
    ZEMd  = np.zeros(nSamples+1)  # Zero Effort Miss distance
    Thoff = np.zeros(nSamples+1)  # Target horiz. offset angle sines
    Tvoff = np.zeros(nSamples+1)  # Target vert. offset angle sines
    

def collectData(i, t, S):
    
    global LastPorY, LastWsgn, LastAsgn
   
    # Get current target and missile states.
    Vt, Pt, Vm, Pm = getS(S)
    
    # Get current missile azimuth and elevation.
    maz, mel = az_el_of_V(Vm)
    
    # Get current target offset horiz and vert angles in
    # missile body frame.
    
    Mbi  = Mrot(maz*RPD, mel*RPD, 0.0)
    Rlos = Prel(Pt, Pm)
    Rtmb = np.matmul(Mbi, Rlos)
    horiz, vert = az_el_of_V(Rtmb)
    
    # Get current target and missile state derivatives.
    dS = dotS(nSvar, S)
    dVt, dPt, dVm, dPm = getSdot(dS)
            
    # Calculate current time-to-go, LOS rate, closing velocity,
    # closing distance, missile acceleration in g's, and ZEM.

    Ulos  = Uvec(Rlos)
    tgo   = timeToGo(Rlos, Vt, Vm)
    Wlosi = Wlos(Vt, Vm, Rlos, Ulos)
    Wlosb = np.matmul(Mbi, Wlosi)
    Amb   = np.matmul(Mbi,dVm)
    if abs(Wlosb[2]) >= abs(Wlosb[1]):
        # predominantly yaw maneuver
        Uzb  = np.matmul(Mbi,  Uzi)  # positive yaw axis
        Uyb  = np.matmul(Mbi,  Uyi)  # positive accel axis
        wsgn = np.sign(np.dot(Wlosb, Uzb))
        asgn = np.sign(np.dot(Amb,   Uyb))
        PorY = 'Y'
    else:
        # predominantly pitch maneuver
        Uyb  = np.matmul(Mbi, -Uyi)  # positive pitch axis
        Uzb  = np.matmul(Mbi,  Uzi)  # positive accel axis
        wsgn = np.sign(np.dot(Wlosb, Uyb))
        asgn = np.sign(np.dot(Amb,   Uzb))
        PorY = 'P'
    losr  = wsgn*la.norm(Wlosb)
    vcls  = la.norm(Vclose(Vt, Vm, Ulos, True))
    dcls  = Dclose(S)
    acmg  = asgn*la.norm(dVm)/g
    zemd  = la.norm(ZEMn(Rlos, Vrel(Vt,Vm), tgo))
    
    Vc, Tgo = calcVcTgo(Pt, Vt, Pm, Vm)
    
    np.testing.assert_almost_equal(vcls, Vc, 3)
    np.testing.assert_almost_equal(tgo, Tgo, 4)
             
    # Display current missile and target states and
    # intercept status; collect data for plotting
    # or printing.
        
    if i == 0 or PRINT_DATA:
        
        print("\nt:  %9.5f" % t)
        print("Target position:   [%9.2f, %9.2f, %9.2f]" % 
              (Pt[0], Pt[1], Pt[2]))
        print("Target velocity:   [%9.3f, %9.3f, %9.3f]" % 
              (Vt[0], Vt[1], Vt[2]))
        print("Missile position:  [%9.2f, %9.2f, %9.2f]" % 
              (Pm[0], Pm[1], Pm[2]))
        print("Missile velocity:  [%9.3f, %9.3f, %9.3f]" % 
              (Vm[0], Vm[1], Vm[2]))
        print("Target velocity magnitude:   %9.3f" % (la.norm(Vt)))
        print("Missile velocity magnitude:  %9.3f" % (la.norm(Vm)))
        print("Missile inertial (az, el):    (%8.3f, %8.3f) degrees" % 
               (maz, mel))
        print("Target offset (horiz, vert):  (%8.3f, %8.3f) degrees" % 
               (horiz, vert))
        print("time to go:  %9.5f sec" % tgo)
        print("LOS rate:    %9.4f deg/sec" % (losr*DPR))
        print("closing velocity:  %9.4f meters/sec" % vcls)
        print("closing distance:  %9.4f meters" % dcls)
        print("missile accel.:    %9.4f g's" % acmg)
        print("ZEM distance:      %9.4f meters" % zemd)
        
    if (PorY != LastPorY) \
        or ((PorY == LastPorY) and ((wsgn != LastWsgn) or (asgn != LastAsgn))):
        print("\nt=%9.5f  %s  sgn_Wlos= %4.1f  sgn_Am= %4.1f" % (t, PorY, wsgn, asgn))
        if PorY != LastPorY : LastPorY = PorY
        if wsgn != LastWsgn : LastWsgn = wsgn
        if asgn != LastAsgn : LastAsgn = asgn
        
    if PLOT_DATA or PRINT_TXYZ:
        
        Time[i]  = t
        Ptx[i]   = Pt[0]
        Pty[i]   = Pt[1]
        Ptz[i]   = Pt[2]
        Vtx[i]   = Vt[0]
        Vty[i]   = Vt[1]
        Vtz[i]   = Vt[2]
        Atx[i]   = dVt[0]
        Aty[i]   = dVt[1]
        Atz[i]   = dVt[2]
        Pmx[i]   = Pm[0]
        Pmy[i]   = Pm[1]
        Pmz[i]   = Pm[2]
        Vmx[i]   = Vm[0]
        Vmy[i]   = Vm[1]
        Vmz[i]   = Vm[2]
        Amx[i]   = dVm[0]
        Amy[i]   = dVm[1]
        Amz[i]   = dVm[2]
        Dcls[i]  = dcls
        LOSd[i]  = losr*DPR
        VELc[i]  = vcls
        Acmg[i]  = acmg
        Velm[i]  = la.norm(Vm)
        ZEMd[i]  = zemd
        Thoff[i] = sin(horiz*RPD)
        Tvoff[i] = sin(vert*RPD)
        
    return


###
### Proportional Navigation main routine
###

if __name__ == "__main__":
    
    # Globals for saving last predominant missile rotation, Pitch
    # or Yaw, and last LOS rate and missile acceleration signs to
    # explain sign changes in LOS rate and missile acceleration 
    # profile plots.
    
    global LastPorY  # pitch or yaw:  'P' or 'Y'
    global LastWsgn  # LOS rate sign: -1.0 or 1.0
    global LastAsgn  # missile acceleration sign: -1.0 or 1.0
    
    LastPorY = None
    LastWsgn = None
    LastAsgn = None
    
    # Estimate time-to-intercept and target position
    # at intercept using Law of Cosines:
    #
    #   C**2 = A**2 + B**2 - 2*A*B*cos(alpha)
    #
    # where:
    #
    #   A     = |Pt0 - Pm0|
    #   B     = |Vm0|*tint
    #   C     = |Vt0|*tint
    #   alpha = leadAngle
    #
    # Then solve the quadratic equation for tint.
    
    magVt = la.norm(Vt0)
    dtm   = la.norm(Prel(Pt0,Pm0))
    alpha = leadAngle(Pt0, Vt0, Pm0, magVm)
    
    c2 = magVm**2 - magVt**2
    c1 = -2*magVm*dtm*cos(alpha)
    c0 = dtm**2
    p  = np.poly1d([c2,c1,c0])
    
    tint  = np.amin(np.abs(p.roots))
    EstPt = Pt0 + Vt0*tint
    
    if tint > tStop:
        print("Error:  estimated time-to-intercept=%8.4f > T_STOP=%8.4f" %
              (tint, tStop))
        sys.exit()
    
    # Calculate azimuth and elevation of estimated target position
    # with respect to missile initial position at time of intercept.
    
    if maz is None:
        maz = atan((EstPt[1]-Pm0[1])/la.norm([EstPt[0]-Pm0[0]]))*DPR
    if mel is None:
        mel = atan((EstPt[2]-Pm0[2])/la.norm([EstPt[0]-Pm0[0], EstPt[1]-Pm0[1]]))*DPR

    Vm0 = setVm(magVm, maz*RPD, mel*RPD)

    
    print("Applying %s proportional navigation guidance law." % (PN_LAWS[PNAV]))
    
    if PRINT_DATA:
        print("Target initial position:  [%9.2f, %9.2f, %9.2f]" % 
              (Pt0[0], Pt0[1], Pt0[2]))
        print("Target initial velocity:  [%9.3f, %9.3f, %9.3f]" % 
              (Vt0[0], Vt0[1], Vt0[2]))
        print("Missile initial position: [%9.2f, %9.2f, %9.2f]" % 
              (Pm0[0], Pm0[1], Pm0[2]))
        print("Missile velocity magnitude:   %9.3f" % magVm)
        print("Estimated time-to-intercept:  %9.5f" % tint)
        print("Estimated target position:  [%9.2f, %9.2f, %9.2f]" % 
              (EstPt[0], EstPt[1], EstPt[2]))
        print("Missile (az, el) lead angles:  (%8.3f, %8.3f) degrees" % 
              (maz, mel))
        print("Missile initial velocity:   [%9.3f, %9.3f, %9.3f]" % 
              (Vm0[0], Vm0[1], Vm0[2]))
    
    # Initialize state vector.
    
    S = setS(S, Vt0, Pt0, Vm0, Pm0)
    
    # Instantiate Runge-Kutta 4th order ode solver, initialize
    # using current state with state time set to zero.
    
    rk4 = RK4_Solver(tStep, nSvar)
    
    S[0] = 0.0
    rk4.init(S)
    
    ## Loop until termination event or simulation stop condition reached.
    
    t = S[0]
    i = 0
    
    while (not Stop(S)) and (i < nSamples):
        
        if i == 0 or PRINT_DATA or PLOT_DATA or PRINT_TXYZ:            
            collectData(i, t, S)
        
        # Perform RK4 integration for animation frame.
        for n in range(0, N_STEP):
            if not Stop(S):
                S = rk4.step(S, dotS)
            
        # Update simulation time and data samples index.
        t = S[0]
        i = i + 1
        
    ## Exited simulation loop.
    
    # Using states and derivatives prior to
    # stop time, reduce step size and integrate
    # till final stop condition.
    
    num_trys = 0  # account for non-covergent intercepts
    max_trys = 4
    
    while (tStep > 0.00005) and (num_trys < max_trys):
        
        print("\nt=%9.5f  tStep=%9.5f  num_trys=%d  Dclose=%f" % 
              (S[0], tStep, num_trys, Dclose(S)))
        
        S = rk4.get_Sprev()
        
        last_tStep = tStep
        tStep      = delT(S)
         
        rk4 = RK4_Solver(tStep, nSvar)
    
        rk4.init(S)
        
        while (not Stop(S)):
            S = rk4.step(S, dotS)
            
        if tStep == last_tStep:
            num_trys =num_trys + 1
    
    S = rk4.get_Sprev()
    t = S[0]
    
    if num_trys == max_trys:
        
        if t < tStop:
            print("\n*** Non-convergent Intercept ***")
        else:
            print("\n*** Insufficient Simulation Time ***")
        
        INTERCEPT = False
        
    else:
    
        INTERCEPT = True
        
    # Display last missile and target states and
    # intercept status; collect data for plotting
    # or printing.
    
    PRINT_DATA = True
    
    collectData(i, t, S)
    
    istop = i
    iend  = istop + 1
    
    if (num_trys == max_trys) or (Dcls[istop] > MinMissDist):
        INTERCEPT = False
        print("\n*** Missile missed target.")
    else:
        INTERCEPT = True
        print("\n*** Missile intercepted target.")
    
    if PLOT_DATA or PRINT_TXYZ :
        
        # Ensure Time data array contains nSamples of simulation time steps
        # in order to prevent plotted lines from wrapping back to time 0.0
        # if the simulation loop was terminated before all nSamples of data
        # were collected.
        
        while i < nSamples:
            Time[i] = t
            t = t + N_TIME
            i = i + 1
     
    if PLOT_DATA:
        
        # Create and show the plot figures.
        
        figures = []
        
        ## Figure 1 - Closing distance at tStop.
        figures.append(plt.figure(1, figsize=(6,3), dpi=80))
        text = "Closing distance ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('Distance (meters)')
        plt.xlim([Time[istop-3], Time[istop]])
        plt.ylim([-1.0, Dcls[istop-3]])
        plt.grid()
        plt.plot(Time[istop-3:iend], Dcls[istop-3:iend], 'o:k')
        plt.plot(np.array([Time[istop-3], Time[istop]]), 
                 np.array([MinMissDist, MinMissDist]), '-c')
        if INTERCEPT:
            plt.plot(Time[istop], Dcls[istop], 'X:m')
            plt.legend(['Distance','MinMissDist','Intercept'])
        else:
            plt.plot(Time[istop], Dcls[istop], 'o:c')
            plt.legend(['Distance','MinMissDist','Missed'])
        
        ## Figure 2 - XY plan view of intercept geometry at tStop.
        figures.append(plt.figure(2, figsize=(6,6), dpi=80))
        text = "XY plan view of intercept ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('X (meters)')
        plt.ylabel('y (meters)')
        if MSL == SAM:
            plt.xlim([(Pmx[istop]+Ptx[istop])/2-10.0,
                      (Pmx[istop]+Ptx[istop])/2+5.0])
            plt.ylim([(Pmy[istop]+Pty[istop])/2-10.0,
                      (Pmy[istop]+Pty[istop])/2+5.0])
        else:
            plt.xlim([(Pmx[istop]+Ptx[istop])/2-20.0,
                      (Pmx[istop]+Ptx[istop])/2+10.0])
            plt.ylim([(Pmy[istop]+Pty[istop])/2-20.0,
                      (Pmy[istop]+Pty[istop])/2+10.0])
        plt.grid()
        plt.plot(Ptx[istop-3:iend], Pty[istop-3:iend], 's:r',
                 Pmx[istop-3:iend], Pmy[istop-3:iend], 'x:b',
                 np.array([Pmx[istop-3],Ptx[istop-3]]),
                 np.array([Pmy[istop-3],Pty[istop-3]]), '.:k',
                 np.array([Pmx[istop-2],Ptx[istop-2]]),
                 np.array([Pmy[istop-2],Pty[istop-2]]), '.:k',
                 np.array([Pmx[istop-1],Ptx[istop-1]]),
                 np.array([Pmy[istop-1],Pty[istop-1]]), '.:k')
        if INTERCEPT:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmy[istop],Pty[istop]]), '.:m')
            plt.legend(('Target','Missile','LOS',' ',' ','Intercept'), loc='upper left')
        else:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmy[istop],Pty[istop]]), '.:c')
            plt.legend(('Target','Missile','LOS',' ', ' ','Missed'), loc='upper left')
        
        ## Figure 3 - XZ profile view of intercept geometry at tStop.
        figures.append(plt.figure(3, figsize=(6,3), dpi=80))
        text = "XZ profile view of intercept ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('X (meters)')
        plt.ylabel('Z (meters)')
        if MSL == SAM:
            plt.xlim([(Pmx[istop]+Ptx[istop])/2-10.0,
                      (Pmx[istop]+Ptx[istop])/2+5.0])
            plt.ylim([(Pmz[istop]+Ptz[istop])/2-5.0,
                      (Pmz[istop]+Ptz[istop])/2+2.5])
        else:
            plt.xlim([(Pmx[istop]+Ptx[istop])/2-20.0,
                      (Pmx[istop]+Ptx[istop])/2+10.0])
            plt.ylim([(Pmz[istop]+Ptz[istop])/2-10.0,
                      (Pmz[istop]+Ptz[istop])/2+5.0])
        plt.grid()
        plt.plot(Ptx[istop-3:iend], Ptz[istop-3:iend], 's:r',
                 Pmx[istop-3:iend], Pmz[istop-3:iend], 'x:b',
                 np.array([Pmx[istop-3],Ptx[istop-3]]),
                 np.array([Pmz[istop-3],Ptz[istop-3]]), '.:k',
                 np.array([Pmx[istop-2],Ptx[istop-2]]),
                 np.array([Pmz[istop-2],Ptz[istop-2]]), '.:k',
                 np.array([Pmx[istop-1],Ptx[istop-1]]),
                 np.array([Pmz[istop-1],Ptz[istop-1]]), '.:k')
        if INTERCEPT:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmz[istop],Ptz[istop]]), '.:m')
            plt.legend(('Target','Missile','LOS',' ',' ','Intercept'), loc='upper left')
        else:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmz[istop],Ptz[istop]]), '.:c')
            plt.legend(('Target','Missile','LOS',' ',' ','Missed'), loc='upper left')
        
        ## Figure 4 - XY plan view of missile/target engagement.
        figures.append(plt.figure(4, figsize=(6,6), dpi=80))
        text = "XY plan view of missile/target engagement ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('X (meters)')
        plt.ylabel('Y (meters)')
        plt.xlim([Pm0[0], Pt0[0]+500.0])
        if MSL == SAM:
            plt.ylim([Pm0[1], Pm0[1]+(Pt0[0]+500-Pm0[0])])
        else:
            plt.ylim([Pm0[1]-Pm0[1]/2.0, Pm0[1]+Pm0[1]/2])
        plt.grid()
        plt.plot(Ptx[0:iend], Pty[0:iend], '-r',
                 Pmx[0:iend], Pmy[0:iend], '-b')
        if INTERCEPT:
            plt.plot(Pmx[istop:iend], Pmy[istop:iend], 'xm')
            plt.legend(('Target','Missile','Intercept'), loc='upper left')
        else:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmy[istop],Pty[istop]]), 'oc')
            plt.legend(('Target','Missile','Missed'), loc='upper left')
        
        ## Figure 5 - XZ profile view of missile/target engagement.
        figures.append(plt.figure(5, figsize=(6,3), dpi=80))
        text = "XZ profile view of missile/target engagement ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('X (meters)')
        plt.ylabel('Z (meters)')
        plt.xlim([Pm0[0], Pt0[0]+500.0])
        if MSL == SAM:
            plt.ylim([0.0, (Pt0[0]+500-Pm0[0])/2])
        else:
            plt.ylim([floor(min(Pmz[0:istop])/1000.0)*1000.0 - 1000.0,
                      floor(min(Pmz[0:istop])/1000.0)*1000.0 + 
                      (Pt0[0]-Pm0[0])/2 - 1000.0])
        plt.grid()
        plt.plot(Ptx[0:iend], Ptz[0:iend], '-r',
                 Pmx[0:iend], Pmz[0:iend], '-b')
        if INTERCEPT:
            plt.plot(Pmx[istop:iend], Pmz[istop:iend], 'xm')
            plt.legend(('Target','Missile','Intercept'), loc='upper left')
        else:
            plt.plot(np.array([Pmx[istop],Ptx[istop]]),
                     np.array([Pmz[istop],Ptz[istop]]), 'oc')
            plt.legend(('Target','Missile','Missed'), loc='upper left')
        
        ## Figure 6 - Missile velocity magnitude vs time of flight
        ##            up to point just before flight termination.
        ##            Note; This is a constant velocity missile.
        figures.append(plt.figure(6, figsize=(6,3), dpi=80))
        text = "Missile velocity magnitude profile ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('Velocity (meters/sec)')
        plt.xlim([0.0, T_STOP])
        plt.ylim([Velm[0]-5, ceil(max(Velm[0:istop]))+5])
        plt.grid()
        plt.plot(Time[0:istop], Velm[0:istop], '-k')
        
        ## Figure 7 - Missile acceleration vs time of flight up
        ##            to point just before flight termination.
        figures.append(plt.figure(7, figsize=(6,3), dpi=80))
        text = "Missile acceleration profile ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('Acceleration (g)')
        plt.xlim([0.0, T_STOP])
        plt.ylim([floor((min(Acmg[0:istop])-Gmmax[MSL]/2)*10.0)/10.0,
                   ceil((max(Acmg[0:istop])+Gmmax[MSL]/2)*10.0)/10.0])
        plt.grid()
        plt.plot(Time[0:istop], Acmg[0:istop], ',-k')
        
        ## Figure 8 - Line-of-Sight rate vs time of flight up
        ##            to point just before flight termination.
        figures.append(plt.figure(8, figsize=(6,3), dpi=80))
        text = "LOS rate profile ({0}, N={1}, At={2})".\
            format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('LOS rate (deg/sec)')
        plt.xlim([0.0, T_STOP])
        plt.ylim([max([floor(min(LOSd[0:istop]))-0.2,-6]), 
                  min([ ceil(max(LOSd[0:istop]))+0.2, 6])])
        plt.grid()
        plt.plot(Time[0:istop], LOSd[0:istop], '-k')
        
        ## Figure 9 - Closing velocity vs time of flight up
        ##            to point just before flight termination.
        figures.append(plt.figure(9, figsize=(6,3), dpi=80))
        text = "Closing velocity profile ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('Velocity (meters/sec)')
        plt.xlim([0.0, T_STOP])
        plt.ylim([floor(min(VELc[0:istop]))-10.0,
                   ceil(max(VELc[0:istop]))+10.0])
        plt.grid()
        plt.plot(Time[0:istop], VELc[0:istop], '-k')
        
        ## Figure 10 - Zero Effort Miss distance vs time of flight
        ##             up to point just before flight termination.
        figures.append(plt.figure(10, figsize=(6,3), dpi=80))
        text = "Zero Effort Miss distance profile ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Time (sec)')
        plt.ylabel('Distance (meters)')
        plt.xlim([0.0, T_STOP])
        plt.ylim([0.0, ceil(max(ZEMd[0:istop]))+10.0])
        plt.grid()
        plt.plot(Time[0:istop], ZEMd[0:istop], '-k')
        
        ## Figure 11 - Target offset sines wrt missile +x axis.
        figures.append(plt.figure(11, figsize=(6,6), dpi=80))
        text = "Target offset sines wrt missile +x axis ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        plt.title(text)
        plt.xlabel('Horizontal Offset (Sine)')
        plt.ylabel('Vertical Offset (Sine)')
        plt.xlim([-0.5,  0.5])  # +/- 30 degrees
        plt.ylim([ 0.5, -0.5])  # +/- 30 degrees
        plt.grid()
        plt.plot(Thoff[0:1], Tvoff[0:1], 'o:r')
        plt.plot(Thoff[1:iend], Tvoff[1:iend], ',:r')
        if INTERCEPT:
            plt.plot(Thoff[istop:iend], Tvoff[istop:iend], 'xm')
            plt.legend(('Target at To','Target Path','Intercept'), loc='upper left')
        else:
            plt.plot(Thoff[istop:iend], Tvoff[istop:iend], 'oc')
            plt.legend(('Target at To','Target Path','Missed'), loc='upper left')
        
        ## Figure 12 - 3D missile/target engagement trajectories plot.
        figures.append(plt.figure(12, figsize=(8,8), dpi=80))
        ax = figures[-1].add_subplot(projection='3d')
        text = "3D Plot of missile/target engagement ({0}, N={1}, At={2})"\
            .format(PN_LAWS[PNAV], int(Nm), int(Nt))
        ax.set_title(text)
        ax.set_xlabel('X (meters)')
        ax.set_ylabel('Y (meters)')
        ax.set_zlabel('Z (meters)')
        if MSL == SAM:
            ax.set_xlim3d([Pm0[0],  Pt0[0]+500.0])
            ax.set_ylim3d([Pm0[1],  Pm0[1]+(Pt0[0]+500-Pm0[0])])
            ax.set_zlim3d([0.0,    (Pt0[0]+500-Pm0[0])/2])
        else:
            ax.set_xlim3d([floor(Pmx[0]/500.0)*500.0,
                           ceil(Ptx[0]/500.0)*500.0])
            ax.set_ylim3d([floor(Pmy[0]/500.0)*500.0 - 
                           (ceil(Ptx[0]/500.0)*500.0 - 
                            floor(Pmx[0]/500.0)*500.0)/2,
                           floor(Pmy[0]/500.0)*500.0 +
                           (ceil(Ptx[0]/500.0)*500.0 - 
                            floor(Pmx[0]/500.0)*500.0)/2])
            ax.set_zlim3d([floor(Pmz[0]/500.0)*500.0 - 
                           (ceil(Ptx[0]/500.0)*500.0 - 
                            floor(Pmx[0]/500.0)*500.0)/2,
                           floor(Pmz[0]/500.0)*500.0 +
                            (ceil(Ptx[0]/500.0)*500.0 - 
                             floor(Pmx[0]/500.0)*500.0)/2])
        ax.grid()
        # target in 3D space
        line1 = Line3D(Ptx[0:iend], Pty[0:iend], Ptz[0:iend],
                       color='r', ls='-', lw=2.0,       
                       marker=' ', mew=2.0, mec='r', mfc='r',
                       label='Target Flight Path')
        # missile in 3D space
        line2 = Line3D(Pmx[0:iend], Pmy[0:iend], Pmz[0:iend],
                       color='b', ls='-', lw=2.0,       
                       marker=' ', mew=2.0, mec='b', mfc='b',
                       label='Missile Trajectory')
        # estimated intercept trajectory in 3D space
        line3 = Line3D(np.array([Pm0[0],EstPt[0]]),
                       np.array([Pm0[1],EstPt[1]]),
                       np.array([Pm0[2],EstPt[2]]),
                       color='g', ls=':', lw=2.0,       
                       marker=' ', mew=1.0, mec='g', mfc='g',
                       label='Estimated Trajectory')
        # estimated intercept point in 3D space
        line4 = Line3D(np.array([EstPt[0],EstPt[0]]),
                       np.array([EstPt[1],EstPt[1]]),
                       np.array([EstPt[2],EstPt[2]]),
                       color='g', ls=' ', lw=1.0,       
                       marker='X', mew=1.0, mec='g', mfc='g',
                       label='Estimated Intercept')
        if INTERCEPT:    
            # intercept point in 3D space
            line5 = Line3D(Pmx[istop:iend], Pmy[istop:iend], Pmz[istop:iend],
                           color='m', ls=' ', lw=1.0,       
                           marker='X', mew=1.0, mec='m', mfc='m',
                           label='Actual Intercept')
        else:
            # missed target and missile points in 3D space
            line5 = Line3D(np.array([Pmx[istop],Ptx[istop]]),
                           np.array([Pmy[istop],Pty[istop]]),
                           np.array([Pmz[istop],Ptz[istop]]),
                           color='c', ls=' ', lw=1.0,       
                           marker='o', mew=1.0, mec='c', mfc='c',
                           label='Missed Intercept')
        ax.add_line(line1)
        ax.add_line(line2)
        ax.add_line(line3)
        ax.add_line(line4)
        ax.add_line(line5)
        ax.legend(loc='upper right')
        ax.view_init(elev=20.0, azim=220.0)
        
        def move_fig(fig):
            """
            Moves given figure plot window based on figure's number.
            """
            fign = fig.number
            x, y = 60*(fign+1), 40*(fign+1)
            backend = mpl.get_backend().upper()
            if backend[0:2] == 'WX':
                fig.canvas.manager.window.SetPosition((x,y))
            elif backend[0:2] == 'TK':
                fig.canvas.manager.window.wm_geometry("+%d+%d" % (x,y)) 
            else:  # QT or GTK
                fig.canvas.manager.window.move(x,y)
        
        try:
            for fig in figures:
                move_fig(fig)
            # Block to keep plots displayed when not running interactively.
            plt.show(block=True)
            plt.close('all')
        except AttributeError:
            # Possibly running "inline".
            plt.show()
        
    if PRINT_TXYZ:
            
        # Set constant missile spin rate. This is purely a 
        # cosmetic effect to provide rendering realism for
        # spin stabilized missiles. Spin is not considered
        # in kinematic equations of motion.

        twopi = 2.0*pi
        fspin = 16.0
        pspin = 1.0/fspin
        wspin = fspin*twopi
        
        # Extract arrays of target and missile position, velocity
        # and acceleration vectors from arrays of saved data. 

        Pte = np.array([[Ptx[0:iend]],
                        [Pty[0:iend]],
                        [Ptz[0:iend]]]).reshape([1,3,iend])
        Vte = np.array([[Vtx[0:iend]],
                        [Vty[0:iend]],
                        [Vtz[0:iend]]]).reshape([1,3,iend])
        Ate = np.array([[Atx[0:iend]],
                        [Aty[0:iend]],
                        [Atz[0:iend]]]).reshape([1,3,iend])
        Pme = np.array([[Pmx[0:iend]],
                        [Pmy[0:iend]],
                        [Pmz[0:iend]]]).reshape([1,3,iend])
        Vme = np.array([[Vmx[0:iend]],
                        [Vmy[0:iend]],
                        [Vmz[0:iend]]]).reshape([1,3,iend])
        Ame = np.array([[Amx[0:iend]],
                        [Amy[0:iend]],
                        [Amz[0:iend]]]).reshape([1,3,iend]) 
    
        #  Calculate missile and target yaw (PSI), pitch (THT) and
        #  roll (PHI) angles.
        #
        #  Note: Vectorized atan2 expression from pitchAngle()
        #        routine used for calculating THTm and THTt.
        
        UVm  = Uvec(Vme[:,0:iend])
        PSIm  = np.arctan2(-UVm[0,1,0:iend], UVm[0,0,0:iend],)*DPR
        THTm  = np.arctan2(UVm[0,2,0:iend], 
                           la.norm([UVm[0,0,0:iend], 
                                    UVm[0,1,0:iend]], axis=0))*DPR
        PHIm  = wspin*Time[0:iend]
        PHImd = np.fmod(PHIm[0:iend], twopi)*DPR
        
        UVt   = Uvec(Vte[:,0:iend])
        PSIt  = np.arctan2(-UVt[0,1,0:iend], UVt[0,0,0:iend])*DPR 
        THTt  = np.arctan2(UVt[0,2,0:iend],
                           la.norm([UVt[0,0,0:iend], 
                                    UVt[0,1,0:iend]], axis=0))*DPR
        PHIt = np.zeros([iend])
        for i in range(0, iend):
            PHIt[i] = bankAngle(Ate[0,:,i],Vte[0,:,i])
        PHItd = np.fmod(PHIt[0:iend], twopi)*DPR
        
        # Assemble run case identifier appended to "./out/TXYZ.OUT."

        sbuff = StringIO()
        sbuff.write("%1d%1d%1d%1d" % (int(MSL), int(PNAV), int(Nm), int(Nt)))
        case = sbuff.getvalue()
        
        # TXYZ output file record formats.

        fmt1 = " %9.4f %9d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f"
        fmt2 = " %9d %9d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f"
        
        # Open TXYZ output file and write trajectory data.
        
        txyz_path = "./out/TXYZ.OUT.{0}".format(case)
        
        try:
            
            f = open(txyz_path, "w")

            for i in range(0, iend):
                print(fmt1 % \
                      (Time[i], 0, Pme[0,0,i], -Pme[0,1,i], -Pme[0,2,i], 
                                   Pte[0,0,i], -Pte[0,1,i], -Pte[0,2,i]),
                      file=f)
                print(fmt2 % \
                      (-9999, -9999, PHImd[i], THTm[i], PSIm[i], 
                                     PHItd[i], THTt[i], PSIt[i]),
                      file=f)
        
            # Append padded time record for rendering and display of
            # last frame by threeD program. Time of istop frame plus 
            # 2*F_TIME, but with missile and target states of istop 
            # frame.
        
            print(fmt1 % \
                  (Time[istop]+2*F_TIME, -1, Pme[0,0,istop], -Pme[0,1,istop], -Pme[0,2,istop], 
                                             Pte[0,0,istop], -Pte[0,1,istop], -Pte[0,2,istop]),
                  file=f)
            print(fmt2 % \
                  (-9999, -9999, PHImd[istop], THTm[istop], PSIm[istop], 
                                 PHItd[istop], THTt[istop], PSIt[istop]),
                  file=f)
            
            print("\nTXYZ.OUT trajectory data written to:  %s" % txyz_path)
                  
            f.close()
            
        except: 
        
            print("\nError:  TXYZ.OUT trajectory data could not be written.")
        