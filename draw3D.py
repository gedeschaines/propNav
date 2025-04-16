# -*- coding: utf-8 -*-

# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name
# pylint: disable=anomalous-backslash-in-string,bad-continuation
# pylint: disable=multiple-statements,redefined-outer-name,global-statement

"""
FILE:  draw3D.py
DATE:  06 DEC 2023
AUTH:  G. E. Deschaines
DESC:  Routines to load, transform and draw polygon data describing
       objects in a 3D world space.
REFS:
     
  [1] This Python script was refactored from draw3D.c available at:
      https://github.com/gedeschaines/threeD/blob/master/src/draw3D.c    
    
NOTE:  Cartesian coordinate frames for world space, field-of-view
       (FOV) viewport (size WxH pixels; aspect ratio (AR) of W/H),
       viewport clipping frustum (pyramid) and drawable pixmap
       are depicted in the following pictograms.
 
                    +X                            +x
          [0,0,0]   /                   [0,0,0]   /  [-W*AR/2,+W*AR/2]
             |     /                       |     /           |
             \--> + ----- +Y               \--> + ----- +y <-/
                  |                             |
                  |                             |   [-H/2,+H/2]
                 +Z                            +z <-----/
  
             World Space                    FOV Viewport
  
  
         [-H/2,+H/2]
             /-> +y  +z                       [0,0]
                  |  /                          + ----- +x [W]
                  | /  [-W/AR/2,+W/AR/2]        |
                  + ----- +x <-/                |
               [0,0,0]                         +y [H]
  
             Clipping Pyramid               Drawable Pixmap

        The x coordinate of points in the FOV viewport are not normalized
        between the FOV frustum near and far planes as in OpenGL or other
        graphics libraries.

Disclaimer:  See DISCLAIMER file.

"""

# Python 3 Module Imports

import sys
import time
from math import floor, sin, cos, tan, asin, atan, atan2, sqrt
from io import StringIO
from copy import copy, deepcopy

# Numpy and Matplotlib Module Imports

try:
    import numpy as np
except ImportError:
    print("* Error: NumPy package required.")
    print("         Suggest installing the SciPy stack.")
    sys.exit()
  
try:
    import matplotlib as mpl
    import matplotlib.collections as mcollections
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath
except ImportError:
    print("* Error: Matplotlib package required.")
    print("         Suggest installing the SciPy stack.")
    sys.exit()


DBG_LVL = 0

# Draw3D Helper Functions (hold over from threeD.c)

def lroundd(x):
    return round(x)

def dmin(x1, x2):
    return min(x1, x2)

def dmax(x1, x2):
    return max(x1, x2)

def lmin(x1, x2):
    return min(x1, x2)

def lmax(x1, x2):
    return max(x1, x2)

# Draw3D Data Structures (using crude implementation of C struct)
        
class Pnt3D:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.X = x
        self.Y = y
        self.Z = z

class PolRec:
    def __init__(self, pt0=Pnt3D(), pt1=Pnt3D(), pt2=Pnt3D()):
        self.Pt0 = pt0
        self.Pt1 = pt1
        self.Pt2 = pt2
        
class Pol3D:
    def __init__(self, k, flg=False, pri=0, typ=0, vis=0, pat=0, \
                 cnt0=Pnt3D(), cnt1=Pnt3D(), nrm0=Pnt3D(), nrm1=Pnt3D()):
        self.Id   = k
        self.Flg  = flg
        self.Pri  = pri
        self.Typ  = typ
        self.Vis  = vis
        self.Pat  = pat
        self.Cnt0 = cnt0
        self.Cnt1 = cnt1
        self.Nrm0 = nrm0
        self.Nrm1 = nrm1

        self.Recs = []
        self.Poly = None

# ThreeD Module Imports

from pqueLib import HeapElement, PriorityQueue 
from clipLib import mxvcnt, polyClip

# Draw3D Constants

fZero = 0.0
fHalf = 0.5
fOne  = 1.0
fTwo  = 2.0
f1K   = 1000.0
rpd   = atan(1.0)/45.0  # radians per degree
dpr   = 1.0/rpd         # degrees per radian

White  = 0
Black  = 1
Red    = 2
Green  = 3
Blue   = 4
Cyan   = 5
Yellow = 6
Brown  = 7
Colors = ['w', 'k', 'r', 'g', 'b', 'c', 'y', 'tab:brown']

# Draw3D Class

class Draw3D:
    
    def __init__(self, txyzFile="TXYZ.OUT", mslType=1, w=600, h=600, 
                       fig=None, ax=None, imgSave=False, imgFPS=50):
        
        # TXYZ file information.
        
        self.txyzFile = txyzFile
        self.mslType = mslType
        self.lfnt = None
        
        # FOV information.

        self.width  = w
        self.height = h
        self.fova   = 90.0      # Field-of-View whole angle (deg)
        self.fovs   = float(h)  # Field-of-View display size (pixels)
        
        # Initialize viewport.
        
        self.xAspect = self.width
        self.yAspect = self.height
        self.ratio   = (self.xAspect*fOne)/(self.yAspect*fOne)
        self.xMax    = lroundd(self.fovs*self.ratio)
        self.yMax    = lroundd(self.fovs)
        self.fovcx   = self.xMax/fTwo
        self.fovcy   = self.yMax/fTwo
        self.fovpt   = Pnt3D()
        
        if DBG_LVL > 0:
            print("Draw3D.__init__:  xAspect,yAspect,ratio = %f %f %f" % \
                  (self.xAspect, self.yAspect, self.ratio))
            print("Draw3D.__init__:  xMax,yMax = %d %d" % \
                  (self.xMax, self.yMax))
            print("Draw3D.__init__:  fovcx,fovcy = %f %f" % \
                  (self.fovcx, self.fovcy))
        
        # Compute viewport FOV focal length and scale factors.
        
        self.zoom = fOne
        self.setFOVsfacs()
        
        # Ground plane grid information.

        self.GridPt1 = [Pnt3D(), Pnt3D(), Pnt3D(), Pnt3D()]
        self.GridPt2 = [Pnt3D(), Pnt3D(), Pnt3D(), Pnt3D()]
        self.GridPathCollections = [None, None]

        # Loaded shape polygon array and vertice buffers.

        self.maxpol = 1024  # Maximum number of polygons
        self.maxpnt =   16  # Maximum number of points in loaded polygon

        self.polcnt  = 0
        self.pollist = [Pol3D(k) for k in range(0, self.maxpol+1)]
        self.pntlist = [Pnt3D() for k in range(0, self.maxpnt+1)]    

        # Priority queue information.
     
        self.MaxElements = 1024
        self.polPQ = PriorityQueue(self.MaxElements)

        self.poltyp_gnd = 0  # ground polygon type
        self.poltyp_tgt = 1  # target polygon type
        self.poltyp_msl = 2  # missile polygon type

        # Main loop processing flags.
        
        self.doneflag = False
        self.quitflag = False
        self.paused = False
        self.align_fov_toward_tgt = False
        self.align_fov_toward_msl = False
        self.align_fov_along_head = True
        self.waitmsec = 0
        
        # Plotting information.
        
        self.fig     = fig
        self.ax      = ax
        self.imgSave = imgSave
        self.imgFPS  = imgFPS
        self.canvas  = fig.canvas
        self.backend = mpl.get_backend().upper()
        self.bckgrnd = self.canvas.copy_from_bbox(self.fig.bbox)
        self.zorder  = 0  # artist drawing order
        
        # Rendered Image Trajectory State Annotations.
        
        self.fr_time = None
        
        self.TextArtists = {'time':None, 'zoom':None,
                            'Xm' :None, 'Ym' :None, 'Zm' :None,
                            'PSm':None, 'THm':None, 'PHm':None,
                            'Xt' :None, 'Yt' :None, 'Zt' :None,
                            'PSt':None, 'THt':None, 'PHt':None}
        
        self.TextFormats = {'time':"Time= %8.4f",
                            'zoom':"Zoom= %8.4f",
                            'Xm'  :"Xm= %10.2f",
                            'Ym'  :"Ym= %10.2f",
                            'Zm'  :"Zm= %10.2f",
                            'PSm' :"PSm= %8.3f",
                            'THm' :"THm= %8.3f",
                            'PHm' :"PHm= %8.3f",
                            'Xt'  :"Xt= %10.2f",
                            'Yt'  :"Yt= %10.2f",
                            'Zt'  :"Zt= %10.2f",
                            'PSt' :"PSt= %8.3f",
                            'THt' :"THt= %8.3f",
                            'PHt' :"PHt= %8.3f"}
        
        self.TextCoords = {'time':(10,12), 'zoom':(10,24),
                            'Xm' :(90,12),  'Ym' :(90,24),  'Zm' :(90,36),
                            'PSm':(170,12), 'THm':(170,24), 'PHm':(170,36),
                            'Xt' :(250,12), 'Yt' :(250,24), 'Zt' :(250,36),
                            'PSt':(330,12), 'THt':(330,24), 'PHt':(330,36)}
        
        
    def setFOVsfacs(self):
        """
        Sets the World to FOV space X, Y, Z coordinate scale factors.
        """
        
        zfovr = fTwo*atan(tan(fHalf*self.fova*rpd)/self.zoom)
        tanfv = sin(zfovr/fTwo)/cos(zfovr/fTwo)
        
        self.fl      = (self.fovs/fTwo)/tanfv
        self.flmin   = 0.1*self.fl
        self.sfacx   = fOne  # Not Used
        self.sfacy   = fOne/tanfv
        self.sfacyAR = self.sfacy/self.ratio
        self.sfacz   = fOne/tanfv
        
        if DBG_LVL > 0:
            print("Draw3D.setFOVsfacs:  tanfv,fl = %f %f" % \
                  (tanfv, self.fl))
            print("Draw3D.setFOVsfacs:  sfacx,sfacy,sfacyAR,sfacz = %f %f %f %f" % \
                  (self.sfacx, self.sfacy, self.sfacyAR, self.sfacz))
   

    def MakePol(self, pntcnt, thePri, theTyp, theVis, thePat, offset):
        """
        Makes list of polygon data records from array of polygon data.
        """

        # Increment polygon counter.

        self.polcnt = self.polcnt + 1

        if DBG_LVL > 0:
            print("MakePol:  polcnt = %d" % (self.polcnt))

        # Initialize polygon list entry.

        self.pollist[self.polcnt].Flg = False
        self.pollist[self.polcnt].Pri = thePri*100000000
        self.pollist[self.polcnt].Pat = thePat
        self.pollist[self.polcnt].Typ = theTyp
        self.pollist[self.polcnt].Vis = theVis

        # Initialize first polygon point record.

        newRec = PolRec(pt0=Pnt3D(x=self.pntlist[1].X + offset.X,
                                  y=self.pntlist[1].Y + offset.Y,
                                  z=self.pntlist[1].Z + offset.Z),
                        pt1=Pnt3D(x=self.pntlist[1].X + offset.X,
                                  y=self.pntlist[1].Y + offset.Y,
                                  z=self.pntlist[1].Z + offset.Z),
                        pt2=Pnt3D(x=fZero, y=fZero, z=fZero))

        if DBG_LVL > 2:
            print("MakePol:  vertice # %d =  %f  %f  %f" % \
                 (1, newRec.Pt0.X, newRec.Pt0.Y, newRec.Pt0.Z))
                
        # Add new polygon point record to polygon.
        
        self.pollist[self.polcnt].Recs.append(deepcopy(newRec))
        
        # Initialize polygon centroid summation variables.

        sumP = fOne
        sumX = deepcopy(newRec.Pt0.X)
        sumY = deepcopy(newRec.Pt0.Y)
        sumZ = deepcopy(newRec.Pt0.Z)

        # Initialize polygon point record for each additional point.

        i = 2
        while i <= pntcnt:
            
            # Initialize new polygon point record.

            newRec = PolRec(pt0=Pnt3D(x=self.pntlist[i].X + offset.X,
                                      y=self.pntlist[i].Y + offset.Y,
                                      z=self.pntlist[i].Z + offset.Z),
                            pt1=Pnt3D(x=self.pntlist[i].X + offset.X,
                                      y=self.pntlist[i].Y + offset.Y,
                                      z=self.pntlist[i].Z + offset.Z),
                            pt2=Pnt3D(x=fZero, y=fZero, z=fZero))
            
            if DBG_LVL > 2:
                print("MakePol:  vertice # %d =  %f  %f  %f" % \
                      (i, newRec.Pt0.X, newRec.Pt0.Y, newRec.Pt0.Z))
            
            # Add new polygon point record to polygon. 
            self.pollist[self.polcnt].Recs.append(deepcopy(newRec))
            
            # Increment polygon centroid summation variables.
            sumP = sumP + fOne
            sumX = sumX + deepcopy(newRec.Pt0.X)
            sumY = sumY + deepcopy(newRec.Pt0.Y)
            sumZ = sumZ + deepcopy(newRec.Pt0.Z)
            
            # Increment polygon point counter.
            i = i + 1

        # Calculate polygon centroid.

        self.pollist[self.polcnt].Cnt0 = Pnt3D(x=sumX/sumP,
                                               y=sumY/sumP,
                                               z=sumZ/sumP)
        self.pollist[self.polcnt].Cnt1 = Pnt3D(x=sumX/sumP,
                                               y=sumY/sumP,
                                               z=sumZ/sumP)
        
        if DBG_LVL > 2:
            print("MakePol:  centroid 0 =  %f  %f  %f" % \
                  (self.pollist[self.polcnt].Cnt0.X,
                   self.pollist[self.polcnt].Cnt0.Y,
                   self.pollist[self.polcnt].Cnt0.Z))
            print("MakePol:  centroid 1 =  %f  %f  %f" % \
                  (self.pollist[self.polcnt].Cnt1.X,
                   self.pollist[self.polcnt].Cnt1.Y,
                   self.pollist[self.polcnt].Cnt1.Z))

        # Calculate polygon normal assuming traversal from point 0 
        # to point 1 is in a counterclockwise direction.
        
        V0x = self.pollist[self.polcnt].Recs[0].Pt0.X - self.pollist[self.polcnt].Cnt0.X
        V0y = self.pollist[self.polcnt].Recs[0].Pt0.Y - self.pollist[self.polcnt].Cnt0.Y
        V0z = self.pollist[self.polcnt].Recs[0].Pt0.Z - self.pollist[self.polcnt].Cnt0.Z
        V1x = self.pollist[self.polcnt].Recs[1].Pt0.X - self.pollist[self.polcnt].Cnt0.X
        V1y = self.pollist[self.polcnt].Recs[1].Pt0.Y - self.pollist[self.polcnt].Cnt0.Y
        V1z = self.pollist[self.polcnt].Recs[1].Pt0.Z - self.pollist[self.polcnt].Cnt0.Z
        
        V0  = np.array([V0x, V0y, V0z])
        V1  = np.array([V1x, V1y, V1z])
        V01 = np.cross(V0, V1)
        nrm = V01 / np.linalg.norm(V01)
        
        self.pollist[self.polcnt].Nrm0 = Pnt3D(x=nrm[0],
                                               y=nrm[1],
                                               z=nrm[2])
        self.pollist[self.polcnt].Nrm1 = Pnt3D(x=nrm[0],
                                               y=nrm[1],
                                               z=nrm[2])
        
        if DBG_LVL > 2:
            print("MakePol:  normal 0 =  %f  %f  %f" % \
                  (self.pollist[self.polcnt].Nrm0.X,
                   self.pollist[self.polcnt].Nrm0.Y,
                   self.pollist[self.polcnt].Nrm0.Z))
            print("MakePol:  normal 1 =  %f  %f  %f" % \
                  (self.pollist[self.polcnt].Nrm1.X,
                   self.pollist[self.polcnt].Nrm1.Y,
                   self.pollist[self.polcnt].Nrm1.Z))
                
    def MakeMatrix(self, p, t, r): 
        """
        Computes world space to view space rotation transformation
        matrix elements.
        
        RHS yaw, pitch, roll as p, t, r respectively in radians.
        """
        
        cosp = cos(p)
        sinp = sin(p)
        cost = cos(t)
        sint = sin(t)
        cosr = cos(r)
        sinr = sin(r)
        
        kctcp   = cost*cosp
        kctsp   = cost*sinp
        kstcp   = sint*cosp
        kstsp   = sint*sinp
        kcrsp   = cosr*sinp
        kcrcp   = cosr*cosp
        ksrsp   = sinr*sinp
        ksrcp   = sinr*cosp
        ksrct   = sinr*cost
        kcrct   = cosr*cost
        ksrstcp = sinr*kstcp
        ksrstsp = sinr*kstsp
        kcrstcp = cosr*kstcp
        kcrstsp = cosr*kstsp
        
        self.dcx1 =  kctcp
        self.dcy1 =  kctsp
        self.dcz1 = -sint
        self.dcx2 =  ksrstcp - kcrsp
        self.dcy2 =  ksrstsp + kcrcp
        self.dcz2 =  ksrct
        self.dcx3 =  kcrstcp + ksrsp
        self.dcy3 =  kcrstsp - ksrcp
        self.dcz3 =  kcrct


    def XfrmGrid(self):
        """
        Transforms grid world space coordinates to viewport coordinates.
        """
        
        if self.mslType == 1:
            # SAM engagement scenario case
            self.GridPt1[0].X =  2000.0
            self.GridPt1[0].Y = -2000.0
            self.GridPt1[0].Z =     0.0
            self.GridPt1[1].X =  2000.0
            self.GridPt1[1].Y =  2000.0
            self.GridPt1[1].Z =     0.0
            self.GridPt1[2].X = -2000.0
            self.GridPt1[2].Y =  2000.0
            self.GridPt1[2].Z =     0.0
            self.GridPt1[3].X = -2000.0
            self.GridPt1[3].Y = -2000.0
            self.GridPt1[3].Z =     0.0
        else:
            # AAM engagement scenario case
            self.GridPt1[0].X =  20000.0
            self.GridPt1[0].Y = -16000.0
            self.GridPt1[0].Z =      0.0
            self.GridPt1[1].X =  20000.0
            self.GridPt1[1].Y =   4000.0
            self.GridPt1[1].Z =      0.0
            self.GridPt1[2].X =      0.0
            self.GridPt1[2].Y =   4000.0
            self.GridPt1[2].Z =      0.0
            self.GridPt1[3].X =      0.0
            self.GridPt1[3].Y = -16000.0
            self.GridPt1[3].Z =      0.0
        
        for k in range(0,4):
            # Translate coordinates into viewport FOV reference frame.
            xd = self.GridPt1[k].X - self.fovpt.X
            yd = self.GridPt1[k].Y - self.fovpt.Y
            zd = self.GridPt1[k].Z - self.fovpt.Z
            # Rotate coordinates into viewport reference frame.
            xs = self.dcx1*xd + self.dcy1*yd + self.dcz1*zd
            ys = self.dcx2*xd + self.dcy2*yd + self.dcz2*zd
            zs = self.dcx3*xd + self.dcy3*yd + self.dcz3*zd
            # Apply FOV scaling.
            ys = ys*self.sfacyAR  # account for square clipping frustum base of fovs pixels
            zs = zs*self.sfacz
            # Save scaled viewport coordinates.
            self.GridPt2[k].X = xs
            self.GridPt2[k].Y = ys
            self.GridPt2[k].Z = zs


    def XfrmPoly(self, iPol):
        """
        Transforms polygon world space coordinates to viewport coordinates.
        """
    
        if DBG_LVL > 3:
            print("XfrmPoly:  Polygon # %d" % (iPol))
        
        # Check if polygon surface is visible.
        
        if self.pollist[iPol].Vis == 2:
            nrm1  = self.pollist[iPol].Nrm1
            nrm2x = self.dcx1*nrm1.X + self.dcy1*nrm1.Y + self.dcz1*nrm1.Z
            nrm2y = self.dcx2*nrm1.X + self.dcy2*nrm1.Y + self.dcz2*nrm1.Z
            nrm2z = self.dcx3*nrm1.X + self.dcy3*nrm1.Y + self.dcz3*nrm1.Z
            dotp  = nrm2x*self.dcx1 + nrm2y*self.dcx2 + nrm2z*self.dcx3
            if dotp > 0.0:
                self.pollist[iPol].Flg = False
                return
            
        # Compute screen coordinates for polygon vertice.  
        
        inflag = False
        pcode  = self.pollist[iPol].Pri
        
        for ipnt in range(0, len(self.pollist[iPol].Recs)):
            
            aPolRec = self.pollist[iPol].Recs[ipnt]
            
            if DBG_LVL > 4:
                print("    - 1 vertice  %2d:  %f  %f  %f" % \
                      (ipnt, aPolRec.Pt1.X, aPolRec.Pt1.Y, aPolRec.Pt1.Z))
                
            # Translate coordinates into viewport FOV reference frame.
            xd = aPolRec.Pt1.X - self.fovpt.X
            yd = aPolRec.Pt1.Y - self.fovpt.Y
            zd = aPolRec.Pt1.Z - self.fovpt.Z
            # Rotate coordinates into viewport reference frame.
            xs = self.dcx1*xd + self.dcy1*yd + self.dcz1*zd
            ys = self.dcx2*xd + self.dcy2*yd + self.dcz2*zd
            zs = self.dcx3*xd + self.dcy3*yd + self.dcz3*zd
            # Apply FOV scaling.
            ys = ys*self.sfacyAR  # account for square clipping frustum base of fovs pixels
            zs = zs*self.sfacz
            # Save scaled viewport coordinates.
            aPolRec.Pt2.X = xs
            aPolRec.Pt2.Y = ys
            aPolRec.Pt2.Z = zs
            if DBG_LVL > 4:
                print("    - 2 vertice  %2d:  %f  %f  %f" % \
                      (ipnt, aPolRec.Pt2.X, aPolRec.Pt2.Y, aPolRec.Pt2.Z))         
            self.pollist[iPol].Recs[ipnt] = deepcopy(aPolRec)
            
            # Check if point is in front of view point.
            if xs >= fZero:
                inflag = True
                
        # Enque polygon if at least one vertice is in the viewport.
        
        if inflag and ( not self.polPQ.isFull() ):
            xd    = self.pollist[iPol].Cnt1.X - self.fovpt.X
            yd    = self.pollist[iPol].Cnt1.Y - self.fovpt.Y
            zd    = self.pollist[iPol].Cnt1.Z - self.fovpt.Z
            xs    = self.dcx1*xd + self.dcy1*yd + self.dcz1*zd
            ys    = self.dcx2*xd + self.dcy2*yd + self.dcz2*zd
            zs    = self.dcx3*xd + self.dcy3*yd + self.dcz3*zd
            rs    = sqrt(xs*xs + ys*ys + zs*zs)
            rsmm  = rs*f1K
            irsmm = lroundd(rsmm)
            anElement = HeapElement(pcode + irsmm, iPol)
            if DBG_LVL > 3:
                print("    - element:  %ld  %hd  %hd  %hd  %ld  %f  %f  %f  %f  %ld" % \
                      (anElement.key,
                             anElement.info,
                                 self.pollist[iPol].Typ,
                                      self.pollist[iPol].Vis,
                                          self.pollist[iPol].Pri,
                                               xs, ys, zs, rsmm, irsmm))
            self.polPQ.priorityEnq(deepcopy(anElement))
            if irsmm < 0:
                print("*** XfrmPoly:  lroundd(rsmm) < 0 for iPol=%hd" % (iPol))
                self.quitflag = True
                
        self.pollist[iPol].Flg = inflag


    def MovePoly(self, iPol, px, py, pz):
        """
        Moves polygons in world space.
        """
    
        # Move the polygon centroid.
        xb                        = copy(self.pollist[iPol].Cnt0.X)
        yb                        = copy(self.pollist[iPol].Cnt0.Y)
        zb                        = copy(self.pollist[iPol].Cnt0.Z)
        self.pollist[iPol].Cnt1.X = self.dcx1*xb + self.dcx2*yb + self.dcx3*zb + px
        self.pollist[iPol].Cnt1.Y = self.dcy1*xb + self.dcy2*yb + self.dcy3*zb + py
        self.pollist[iPol].Cnt1.Z = self.dcz1*xb + self.dcz2*yb + self.dcz3*zb + pz
        if DBG_LVL > 3:
            print("MovePoly:  Polygon # %d" % (iPol))
            print("    - centroid :  %f  %f  %f" % \
                                (self.pollist[iPol].Cnt1.X,
                                     self.pollist[iPol].Cnt1.Y,
                                         self.pollist[iPol].Cnt1.Z))
            
        # Move the polygon vertice.
        for ipnt in range(0, len(self.pollist[iPol].Recs)):
            aPolRec       = self.pollist[iPol].Recs[ipnt]
            xb            = copy(aPolRec.Pt0.X)
            yb            = copy(aPolRec.Pt0.Y)
            zb            = copy(aPolRec.Pt0.Z)
            aPolRec.Pt1.X = self.dcx1*xb + self.dcx2*yb + self.dcx3*zb + px
            aPolRec.Pt1.Y = self.dcy1*xb + self.dcy2*yb + self.dcy3*zb + py
            aPolRec.Pt1.Z = self.dcz1*xb + self.dcz2*yb + self.dcz3*zb + pz
            if DBG_LVL > 3:
                print("    - vertice :  %f  %f  %f" % \
                                       (aPolRec.Pt1.X,
                                            aPolRec.Pt1.Y,
                                                aPolRec.Pt1.Z))          
            self.pollist[iPol].Recs[ipnt] = deepcopy(aPolRec)
            
        # Compute moved polygon normal assuming traversal from point 0 
        # to point 1 is in a counterclockwise direction.
        
        V0x = self.pollist[iPol].Recs[0].Pt1.X - self.pollist[iPol].Cnt1.X
        V0y = self.pollist[iPol].Recs[0].Pt1.Y - self.pollist[iPol].Cnt1.Y
        V0z = self.pollist[iPol].Recs[0].Pt1.Z - self.pollist[iPol].Cnt1.Z
        V1x = self.pollist[iPol].Recs[1].Pt1.X - self.pollist[iPol].Cnt1.X
        V1y = self.pollist[iPol].Recs[1].Pt1.Y - self.pollist[iPol].Cnt1.Y
        V1z = self.pollist[iPol].Recs[1].Pt1.Z - self.pollist[iPol].Cnt1.Z
        
        V0  = np.array([V0x, V0y, V0z])
        V1  = np.array([V1x, V1y, V1z])
        V01 = np.cross(V0, V1)
        nrm = V01 / np.linalg.norm(V01)
        
        self.pollist[iPol].Nrm1.X = nrm[0]
        self.pollist[iPol].Nrm1.Y = nrm[1]
        self.pollist[iPol].Nrm1.Z = nrm[2]


    def DrawGrid3D(self, iaxis):
        """
        Draws grid lines clipped to 3D viewing pyramidal frustum.
        """
        
        # Create clipping workspace polygon vertice arrays.
     
        vcnt   = [0, 0, 0, 0, 0, 0, 0, 0]
        vlist  = [[Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)],
                  [Pnt3D() for i in range(0, mxvcnt)]]
        
        # Create clipped path verts and codes arrays.
        
        numLines = 41          # number of grid axis lines
        maxVerts = numLines*2  # max number of grid path vertice
        
        tmpVerts = np.zeros((maxVerts,2))
        tmpCodes = np.zeros(maxVerts)

        if self.GridPathCollections[iaxis-1] is None:
            # Create a temporary path which will be replaced.
            path = mpath.Path([[0.0,        self.height-50.0], 
                               [self.width, self.height-50.0]],
                               [1, 2])
            # Instantiate a PathCollection which will be redrawn.
            self.GridPathCollections[iaxis-1] = \
                 mcollections.PathCollection([path],
                                             animated=False,
                                             visible=False,
                                             ls='-',
                                             lw=1.0,
                                             color=Colors[White])
            self.ax.add_collection(self.GridPathCollections[iaxis-1])
        
        if DBG_LVL > 4:
            print("DrawGrid3D:  Drawing grid axis %d..." % (iaxis))

        if iaxis == 1:
            # Draw grid lines parallel to workd x-axis.
            i10 = 3
            i11 = 2
            i20 = 0
            i21 = 1
        else:
            # Draw grid lines parallel to workd y-axis.
            i10 = 3
            i11 = 0
            i20 = 2
            i21 = 1

        # Calculate incremental distances.

        xd1 = 0.025*(self.GridPt2[i11].X - self.GridPt2[i10].X)
        yd1 = 0.025*(self.GridPt2[i11].Y - self.GridPt2[i10].Y)
        zd1 = 0.025*(self.GridPt2[i11].Z - self.GridPt2[i10].Z)
        xd2 = 0.025*(self.GridPt2[i21].X - self.GridPt2[i20].X)
        yd2 = 0.025*(self.GridPt2[i21].Y - self.GridPt2[i20].Y)
        zd2 = 0.025*(self.GridPt2[i21].Z - self.GridPt2[i20].Z)

        n = 0  # initialize number of path vertice.
        
        for k in range(0, numLines):
            
            pcnt = 1
            
            # Calculate x coordinate of un-clipped line.
            vlist[pcnt][1].X = self.GridPt2[i10].X + k*xd1
            vlist[pcnt][2].X = self.GridPt2[i20].X + k*xd2
             
            if not ((vlist[pcnt][1].X <= self.flmin) and 
                    (vlist[pcnt][2].X <= self.flmin)):
                # Create un-clipped line.
                vlist[pcnt][1].Y = self.GridPt2[i10].Y + k*yd1
                vlist[pcnt][1].Z = self.GridPt2[i10].Z + k*zd1
                vlist[pcnt][2].Y = self.GridPt2[i20].Y + k*yd2
                vlist[pcnt][2].Z = self.GridPt2[i20].Z + k*zd2
                vlist[pcnt][3]   = copy(vlist[pcnt][1])
                vcnt[pcnt]       = 3
                
                if DBG_LVL > 4:
                    for i in range(1,vcnt[pcnt]+1):
                        print("DrawGrid3D:  %d %d %f %f %f" % \
                              (pcnt, i,
                               vlist[pcnt][i].X,
                               vlist[pcnt][i].Y,
                               vlist[pcnt][i].Z))
                
                # Create clipped line.
                pcnt = polyClip(pcnt, vcnt, vlist)
                
                if DBG_LVL > 4:
                    for i in range(1,vcnt[pcnt]+1):
                        print("DrawGrid3D:  %d %d %f %f %f" % \
                              (pcnt, i,
                               vlist[pcnt][i].X,
                               vlist[pcnt][i].Y,
                               vlist[pcnt][i].Z))
                
                # Draw clipped line.
                if vcnt[pcnt] > 2:
                    c = 1  # path code (1:moveto, 2:lineto)
                    for i in range(1,3):
                        xs            = vlist[pcnt][i].X
                        ys            = vlist[pcnt][i].Y/self.sfacyAR
                        zs            = vlist[pcnt][i].Z/self.sfacz
                        sf            = self.fl/xs
                        tmpVerts[n,0] = lroundd(sf*ys) + floor(self.fovcx)
                        tmpVerts[n,1] = lroundd(sf*zs) + floor(self.fovcy)
                        tmpCodes[n]   = c
                        n += 1
                        c = 2
        
        # Reset grid axis PathCollection path vertice and codes, then redraw.        
        self.zorder += fOne
        if n > 1:
            path = mpath.Path(tmpVerts[0:n,:], tmpCodes[0:n])
            self.GridPathCollections[iaxis-1].set_paths([path])
            self.GridPathCollections[iaxis-1].set_visible(True)
            self.GridPathCollections[iaxis-1].set_zorder(self.zorder)
            self.ax.draw_artist(self.GridPathCollections[iaxis-1])
        else:
            # No grid lines to draw, make last grid axis path 
            # collection invisible.
            self.GridPathCollections[iaxis-1].set_visible(False)
        

    def DrawPoly3D(self, iPol):
        """
        Draws polygon clipped to 3D viewing pyramidal frustum.
        """

        # Create clipping workspace polygon vertice arrays.
        
        vcnt  = [0, 0, 0, 0, 0, 0, 0, 0]
        vlist = [[Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)],
                 [Pnt3D() for i in range(0, mxvcnt)]]
        
        # Create clipped polygon XY coordinates array.
        
        xy = np.zeros((mxvcnt,2))

        # Get un-clipped polygon.

        if DBG_LVL > 4:
            print("DrawPoly3D:  Drawing polygon %d..." % (iPol))

        pcnt = 1
        icnt = 0
        for ipnt in range(0, len(self.pollist[iPol].Recs)):
            aPolRec           = self.pollist[iPol].Recs[ipnt]
            icnt              = icnt + 1
            vlist[pcnt][icnt] = aPolRec.Pt2

        icnt              = icnt + 1
        vlist[pcnt][icnt] = vlist[pcnt][1]
        vcnt[pcnt]        = icnt

        if DBG_LVL > 4:
            for i in range(1,vcnt[pcnt]+1):
                print("DrawPoly3D:  %d %d %f %f %f" % \
                      (pcnt, i,
                       vlist[pcnt][i].X,
                       vlist[pcnt][i].Y,
                       vlist[pcnt][i].Z))

        # Create clipped polygon.

        pcnt = polyClip(pcnt, vcnt, vlist)

        if DBG_LVL > 4:
            for i in range(1,vcnt[pcnt]+1):
                print("DrawPoly3D:  %d %d %f %f %f" % \
                      (pcnt, i,
                       vlist[pcnt][i].X,
                       vlist[pcnt][i].Y,
                       vlist[pcnt][i].Z))

        # Draw clipped polygon.

        if vcnt[pcnt] > 3:
            for i in range(1,vcnt[pcnt]+1):
                xs        = vlist[pcnt][i].X
                ys        = vlist[pcnt][i].Y/self.sfacyAR
                zs        = vlist[pcnt][i].Z/self.sfacz
                sf        = self.fl/xs
                xy[i-1,0] = lroundd(sf*ys) + floor(self.fovcx)
                xy[i-1,1] = lroundd(sf*zs) + floor(self.fovcy)
            
            self.zorder += fOne
            if self.pollist[iPol].Vis > 0:
                # Filled polygon.
                if self.pollist[iPol].Poly is None:
                    # Instantiate Patch Polygon which will be redrawn.
                    poly = mpatches.Polygon(xy[0:vcnt[pcnt],:],
                                            closed=True,
                                            animated=False,
                                            visible=True,
                                            fill=True,
                                            lw=1.0,
                                            fc=self.pollist[iPol].Pat,
                                            ec=self.pollist[iPol].Pat,
                                            zorder=self.zorder)
                    self.pollist[iPol].Poly = poly
                    self.ax.add_patch(self.pollist[iPol].Poly)
                else:
                    # Reset the Patch Polygon XY data.
                    self.pollist[iPol].Poly.set_xy(xy[0:vcnt[pcnt],:])
                    self.pollist[iPol].Poly.set(zorder=self.zorder)
                    self.pollist[iPol].Poly.set(visible=True)
            else:
                # Unfilled polygon.
                if self.pollist[iPol].Poly is None:
                    # Instantiate Patch Polygon which will be redrawn.
                    poly = mpatches.Polygon(xy[0:vcnt[pcnt],:],
                                            closed=True,
                                            animated=False,
                                            visible=True,
                                            fill=False,
                                            lw=2.0,
                                            ec=self.pollist[iPol].Pat,
                                            zorder=self.zorder)
                    self.pollist[iPol].Poly = poly
                    self.ax.add_patch(self.pollist[iPol].Poly)
                else:
                    # Reset the Patch Polygon XY data.
                    self.pollist[iPol].Poly.set_xy(xy[0:vcnt[pcnt],:])
                    self.pollist[iPol].Poly.set(zorder=self.zorder)
                    self.pollist[iPol].Poly.set(visible=True)
                    
            # Redraw patch polygon.
            self.ax.draw_artist(self.pollist[iPol].Poly)
      
            
    def LoadPoly(self, lfni, polyfile):
        """
        Read shape polygon data file to make polygon.
        """
        # Read shape model offsets, scaling factor and name record.
        
        line = lfni.readline().strip()
        if line == '':
            print("LoadPoly:  readline error for 1st record in polyfile %s." % \
                  (polyfile))
            return
   
        words = line.split()
        if len(words) < 5:
            print("LoadPoly:  split error for 1st record in polyfile %s." % \
                  (polyfile))
            return
             
        mdloffx = float(words[0])
        mdloffy = float(words[1])
        mdloffz = float(words[2])
        mdlsfc  = float(words[3])
        for i in range(0,len(words)-4):
            if i == 0:
                mdlnam = words[4+i]
            else:
                mdlnam = mdlnam + " " + words[4+i]
    
        # Read shape polygon vertice.

        line = lfni.readline()
    
        while not ((line == '') or (self.polcnt == self.maxpol)):
            
            # Load polygon specification record.
            words = line.strip().split()
            if len(words) >= 7:
                polpnt = int(words[0])
                polpri = int(words[1])
                polcol = int(words[2])
                poltyp = int(words[3])
                polvis = int(words[4])
                polsfc = float(words[5])
                for i in range(0,len(words)-6):
                    if i == 0:
                        polnam = words[6+i]
                    else:
                        polnam = polnam + " " + words[6+i]
                if DBG_LVL > 1:
                    print("LoadPoly:  Loaded specs -  %d  %d  %d  %d  %d %f  %s" % \
                          (polpnt, polpri, polcol, poltyp, polvis, polsfc, polnam))

                # Load vertex points.
                for i in range(1, polpnt+1):
                    line  = lfni.readline().strip()
                    words = line.split()
                    if len(words) == 3:
                        x = float(words[0])
                        y = float(words[1])
                        z = float(words[2])
                        if DBG_LVL > 1:
                            print("LoadPoly:  Loaded vertex -  %f  %f  %f" % \
                                  (x, y, z))
                        self.pntlist[i].X = x*polsfc*mdlsfc
                        self.pntlist[i].Y = y*polsfc*mdlsfc
                        self.pntlist[i].Z = z*polsfc*mdlsfc
                        if DBG_LVL > 1:
                            print("LoadPoly:  Scaled vertex -  %f  %f  %f" % \
                                  (self.pntlist[i].X,
                                   self.pntlist[i].Y,
                                   self.pntlist[i].Z))
                            
                # Load offset point.
                line  = lfni.readline().strip()
                words = line.split()
                if len(words) == 3:
                    x = float(words[0])
                    y = float(words[1])
                    z = float(words[2])
                    if DBG_LVL > 1:
                        print("LoadPoly:  Loaded offset -  %f  %f  %f" % \
                              (x, y, z))
                    offset = Pnt3D(x=x*polsfc*mdlsfc + mdloffx,
                                   y=y*polsfc*mdlsfc + mdloffy,
                                   z=z*polsfc*mdlsfc + mdloffz)
                    if DBG_LVL > 1:
                        print("LoadPoly:  Scaled offset -  %f  %f  %f" % \
                              (offset.X, offset.Y, offset.Z))

                if ( (polcol >= 0) and (polcol < 8) ):
                    self.MakePol(polpnt, polpri, poltyp, polvis, Colors[polcol], deepcopy(offset))
                else:
                    self.MakePol(polpnt, polpri, poltyp, polvis, Colors[Black], deepcopy(offset))

            line = lfni.readline()
    
    
    def ReadPolyData(self):
        """
        Read shape polygon files to load all polygon data.
        """
        
        self.polcnt = 0
        
        if self.mslType == 1:
            lfni = open("./dat/grndpoly1.dat", "rt")
        else:
            lfni = open("./dat/grndpoly2.dat", "rt")
        if lfni is not None:
            if DBG_LVL > 0:
                print("Draw3D:  Loading polygons from file %s" % \
                      ("grndpoly.dat"))
            self.LoadPoly(lfni, "./dat/grndpoly.dat")
            lfni.close()
        
        lfni = open("./dat/fwngpoly.dat", "rt")
        if lfni is not None:
            if DBG_LVL > 0:
                print("Draw3D:  Loading polygons from file %s" % \
                      ("fwngpoly.dat"))
            self.LoadPoly(lfni, "./dat/fwngpoly.dat")
            lfni.close()
        
        if self.mslType == 1:
            lfni = open("./dat/mislpoly1.dat", "rt")
        else:
            lfni = open("./dat/mislpoly2.dat", "rt")
        if lfni is not None:
            if DBG_LVL > 0:
                print("Draw3D:  Loading polygons from file %s" % \
                      ("mislpoly.dat"))
            self.LoadPoly(lfni,"./dat/mislpoly.dat")
            lfni.close()
        

    def AnnotateText(self, text, valnum):
        """
        Creates/updates annotation text artists.
        """
        sbuff = StringIO()
        sbuff.write(self.TextFormats[text] % (valnum))
        valstr =sbuff.getvalue()
        sbuff.close()
        if self.TextArtists[text] is None:
            #self.zorder += fOne
            self.TextArtists[text] = \
                self.ax.annotate(valstr,
                                 self.TextCoords[text],
                                 xycoords="data",
                                 xytext=(0,0),
                                 textcoords=("offset points"),
                                 #backgroundcolor=Colors[Black],
                                 color=Colors[Black],
                                 fontname='monospace',
                                 fontsize='x-small',
                                 ha="left",
                                 va="top",
                                 animated=False,
                                 zorder=self.zorder)
        else:
            self.TextArtists[text].set_text(valstr)
        
    
    def DrawText(self):
        """
        Draws trajectory state values.
        """
        for text in self.TextArtists:
            self.ax.draw_artist(self.TextArtists[text])
        
    
    def onPress(self, event):
        """
        Keyboard key press handler.
        """
        if event.key == 't':
            self.align_fov_toward_tgt = not self.align_fov_toward_tgt
            self.align_fov_toward_msl = False
            self.align_fov_along_head = False
        elif event.key == 'm':
            self.align_fov_toward_msl = not self.align_fov_toward_msl
            self.align_fov_toward_tgt = False
            self.align_fov_along_head = False
        elif event.key == 'h':
            self.align_fov_along_head = True
            self.align_fov_toward_msl = False
            self.align_fov_toward_tgt = False
        elif event.key == 'z':
            self.zoom = fOne
            self.setFOVsfacs()
        elif event.key == 'up':
            self.zoom = self.zoom*1.25
            self.setFOVsfacs()
        elif event.key == 'down':
            self.zoom = self.zoom/1.25
            self.setFOVsfacs()
        elif event.key == '0':
            self.waitmsec = 0
        elif event.key == 'right':
            self.waitmsec -= 10  
            self.waitmsec = lmax(0, self.waitmsec)
        elif event.key == 'left':
            self.waitmsec += 10
            self.waitmsec = lmin(200, self.waitmsec)
        elif event.key == ' ':
            self.paused = not self.paused
        elif event.key == 'x':
            # Exits MainLoop and permits restart.
            self.quitflag = True
            self.doneflag = False
        elif event.key == 'escape':
            # Exits MainLoop and closes figure.
            self.quitflag = True
            self.doneflag = True
        
    
    def MainLoop(self):
        
        self.doneflag = False
        self.quitflag = False
        self.paused = False
        self.align_fov_toward_tgt = False
        self.align_fov_toward_msl = False
        self.align_fov_along_head = True
        self.waitmsec = 0
        
        true_tsec =  0.0
        last_tsec = -1.0/self.imgFPS
        img_count = 0
        img_dtsec = 1.0/self.imgFPS

        # May need to save last missile position.
        
        XM = 0.0 
        YM = 0.0
        ZM = 0.0
        
        # Open trajectory data file.

        if DBG_LVL > 0:
            print("MainLoop:  Opening trajectory file %s" % \
                  (self.txyzFile))
        self.lfnt = open(self.txyzFile, "rt")

        # Main processing loop over trajectory data file.

        line = self.lfnt.readline()
        
        while not ( (line == '') or self.quitflag ):
            
            cpumsec1 = 1000.0*time.perf_counter()
            
            if self.paused: 
                self.canvas.flush_events()
                continue
            
            # Save last "true" missile position (i.e., that
            # read from a previous ktot >= 0 record).
            if true_tsec > 0.0:
                last_XM = XM
                last_YM = YM
                last_ZM = ZM
            
            # Get missile and target position.
            words = line.strip().split()
            if len(words) != 8:
                line = self.lfnt.readline()
                continue           
            tsec = float(words[0])
            ktot = int(words[1]) 
            XM   = float(words[2])
            YM   = float(words[3])
            ZM   = float(words[4])
            XT   = float(words[5])
            YT   = float(words[6])
            ZT   = float(words[7])
            
            if DBG_LVL > 1:
                print("MainLoop:  Time %8.4f, ktot=%4d" % (tsec, ktot))
               
            # Get missile and target orientation.
            line = self.lfnt.readline().strip()
            if line.find("-9999     -9999") < 0:
                # Other 6-DOF simulation trajectory output file.
                words = line.split()
                if len(words) < 6:
                    # Skip this line and process next line.
                    line = self.lfnt.readline()
                    continue
                PHM = float(words[0])
                THM = float(words[1])
                PSM = float(words[2])
                PHT = float(words[3])
                THT = float(words[4])
                PST = float(words[5])
            else:
                # propNav 3-DOF simulation trajectory output file.
                words = line.split()
                if len(words) < 8:
                    # Skip this line and process next line.
                    line = self.lfnt.readline()
                    continue
                PHM  = float(words[2])
                THM  = float(words[3])
                PSM  = float(words[4])
                PHT  = float(words[5])
                THT  = float(words[6])
                PST  = float(words[7])
            
            # skip decoy position and radiance.
            if ktot > 0:
                for itot in range(0, ktot):
                    line = self.lfnt.readline()
            
            # Get target position components. 
            px = XT
            py = YT
            pz = ZT
            
            # Convert target orientation values to radians.
            p = PST*rpd
            t = THT*rpd
            r = PHT*rpd
            if DBG_LVL > 2:
                print("  tgt - px,py,pz,p,t,r = %f %f %f %f %f %f" % \
                      (px, py, pz, p, t, r))
            
            # Compute polygon rotation transformation matrix.
            if DBG_LVL > 2:
                print("MainLoop:  Make target polygon transformation matrix...")
            self.MakeMatrix(p, t, r)
            if DBG_LVL > 2:
                print("           %f  %f  %f" % \
                      (self.dcx1, self.dcy1, self.dcz1))
                print("           %f  %f  %f" % \
                      (self.dcx2, self.dcy2, self.dcz2))
                print("           %f  %f  %f" % \
                      (self.dcx3, self.dcy3, self.dcz3))
            
            # Move target polygons.
            if DBG_LVL > 2:
                print("MainLoop:  Move target polygons...")
            for i in range(1,self.polcnt+1):
                if self.pollist[i].Typ == self.poltyp_tgt:
                    self.MovePoly(i, px, py, pz)
           
            # Get missile position components.
            px = XM
            py = YM
            pz = ZM
            
            # Convert missile orientation values to radians.
            p = PSM*rpd
            t = THM*rpd
            r = PHM*rpd
            if DBG_LVL > 2:
                print("  msl - px,py,pz,p,t,r = %f %f %f %f %f %f" % \
                      (px, py, pz, p, t, r))
                    
            # Compute polygon rotation transformation matrix.
            if DBG_LVL > 2:
                print("MainLoop:  Make missile polygon transformation matrix...")
            self.MakeMatrix(p, t, r)
            if DBG_LVL > 2:
                print("           %f  %f  %f" % \
                      (self.dcx1, self.dcy1, self.dcz1))
                print("           %f  %f  %f" % \
                      (self.dcx2, self.dcy2, self.dcz2))
                print("           %f  %f  %f" % \
                      (self.dcx3, self.dcy3, self.dcz3))
            
            # Move missile polygons.
            if DBG_LVL > 2:
                print("MainLoop:  Move missile polygons...")
            for i in range(1,self.polcnt+1):
                if self.pollist[i].Typ == self.poltyp_msl:
                    self.MovePoly(i, px, py, pz)
            
            # Calculate unit vector from missile to target.
            # NOTE: RHS where +X is forward, +Y is to the
            #       right and +Z is down; -Z is up.
            DXTM = XT - XM  # NOTE: It's highly improbable missile and
            DYTM = YT - YM  #       target positions being identical,
            DZTM = ZT - ZM  #       yielding [DXTM,DYTM,DZTM]=[0,0,0].
            RTM  = sqrt(DXTM*DXTM + DYTM*DYTM + DZTM*DZTM)
            if RTM > 0.0:
                UXTM = DXTM/RTM
                UYTM = DYTM/RTM
                UZTM = DZTM/RTM
            elif ktot > -1:
                # Use last valid missile velocity direction vector.
                DXTM = XM - last_XM
                DYTM = YM - last_YM
                DZTM = ZM - last_ZM
                RTM  = sqrt(DXTM*DXTM + DYTM*DYTM + DZTM*DZTM)
                UXTM = DXTM/RTM
                UYTM = DYTM/RTM
                UZTM = DZTM/RTM
            # Calculate FOV position and orientation.
            if self.align_fov_toward_tgt:
                # Place fovpt near missile; align fov normal axis with
                # unit vector from missile to target.
                self.fovpt.X = XM - fTwo*UXTM
                self.fovpt.Y = YM - fTwo*UYTM
                self.fovpt.Z = dmin(ZM - fTwo*UZTM + fHalf, -0.1)  # keep fovpt above ground
                p = atan2(UYTM, UXTM)  # Yaw    NOTE: Gimbal lock occurs when Pitch is
                t = asin(-UZTM)        # Pitch        +/- 90 deg as both UXTM and UYTM
                r = fZero              # Roll         are zero and Yaw is indeterminate.
            elif self.align_fov_toward_msl:
                # Place fovpt near target; align fov normal axis with
                # unit vector from target to missile.
                self.fovpt.X = XT + 30.0*UXTM
                self.fovpt.Y = YT + 30.0*UYTM
                self.fovpt.Z = ZT + 30.0*UZTM + 15.0
                p = atan2(-UYTM, -UXTM)
                t = asin(UZTM)
                r = fZero
            else:
                # Place fovpt near missile; align fov normal axis with
                # missile heading, but keep in horizontal plane.
                self.fovpt.X = XM - 3.0*cos(p)
                self.fovpt.Y = YM - 3.0*sin(p)
                self.fovpt.Z = dmin(ZM - 1.5, -0.1)  # keep fovpt above ground
                t = fZero
                r = fZero
            
            # Compute FOV rotation transformation matrix.
            if DBG_LVL > 2:
                print("MainLoop:  Make field-of-view rotation matrix...")
            self.MakeMatrix(p, t, r)
            if DBG_LVL > 2:
                print("           %f  %f  %f" % \
                      (self.dcx1, self.dcy1, self.dcz1))
                print("           %f  %f  %f" % \
                      (self.dcx2, self.dcy2, self.dcz2))
                print("           %f  %f  %f" % \
                      (self.dcx3, self.dcy3, self.dcz3))
            
            # Transform ground plane polygon into viewport.
            if DBG_LVL > 2:
                print("MainLoop:  Transform ground plane polygon...")
            self.XfrmPoly(1)
            
            # Transform ground plane grid into viewport.
            if DBG_LVL > 2:
                print("MainLoop:  Transform ground plane grid...")
            self.XfrmGrid()
            
            # Transform object polygons into viewport.
            if DBG_LVL > 2:
                print("MainLoop:  Transform polygons...")
            self.polPQ.clear()
            for i in range(2,self.polcnt+1):
                self.XfrmPoly(i)
            
            # Restore canvas, reset artist drawing zorder, and
            # make all previously drawn polygons invisible.
            self.canvas.restore_region(self.bckgrnd)
            self.zorder = fZero
            for i in range(1, self.polcnt+1):
                if not (self.pollist[i].Poly is None):
                    self.pollist[i].Poly.set(visible=False)
            
            # Draw ground plane polygon.
            if DBG_LVL > 2:
                print("MainLoop:  Draw ground plane polygon...")
            if self.pollist[1].Flg:
                self.DrawPoly3D(1)
            
            # Draw ground plane grid.
            if DBG_LVL > 2:
                print("MainLoop:  Draw ground plane grid...")
            self.DrawGrid3D(1)
            self.DrawGrid3D(2)
            
            # Draw target and missile shape polygons.
            if DBG_LVL > 2:
                print("MainLoop:  Draw target and missile polygons...")
            while not self.polPQ.isEmpty():
                anElement = self.polPQ.priorityDeq()
                if DBG_LVL > 3:
                    print("  %ld  %hd  %hd  %hd  %ld" % \
                          (anElement.key,
                           anElement.info,
                           self.pollist[anElement.info].Typ,
                           self.pollist[anElement.info].Vis,                           
                           self.pollist[anElement.info].Pri))
                self.DrawPoly3D(anElement.info)
            
            # Display time, zoom, missile and target state variables.
            self.zorder = 1000.0
            if ( ktot < 0 ):
                # TXYZ padded time record.
                self.AnnotateText('time', true_tsec)
            else:
                # TXYZ true time record.
                self.AnnotateText('time', tsec)
                true_tsec = tsec
            self.AnnotateText('zoom', self.zoom)
            self.AnnotateText('Xm',  XM)
            self.AnnotateText('Ym',  YM)
            self.AnnotateText('Zm', -ZM)
            self.AnnotateText('PSm', PSM)
            self.AnnotateText('THm', THM)
            self.AnnotateText('PHm', PHM)
            self.AnnotateText('Xt',  XT)
            self.AnnotateText('Yt',  YT)
            self.AnnotateText('Zt', -ZT)
            self.AnnotateText('PSt', PST)
            self.AnnotateText('THt', THT)
            self.AnnotateText('PHt', PHT)
            self.DrawText()
            
            # Update display of rendered image.         
            if self.backend[0:2] == 'WX':
                self.canvas.update()
            elif self.backend[0:2] == 'TK':
                self.canvas.blit(self.ax.bbox)
            elif self.backend[0:2] == 'NB':
                self.canvas.blit(self.ax.bbox)
            elif self.backend == 'MODULE://IPYMPL.BACKEND_NBAGG':
                self.canvas.blit(self.ax.bbox)
            else:  # QT or GTK
                self.canvas.update()
            self.canvas.flush_events()
            
            # Save figure canvas with rendered image to PNG file.
            if self.imgSave:
                if ( (tsec+0.005 - last_tsec) >= img_dtsec ):
                    sbuff = StringIO()
                    sbuff.write("./Ximg/img_%04hd.png" % (img_count))
                    img_fname = sbuff.getvalue()
                    sbuff.close()
                    self.fig.savefig(img_fname, format='png')
                    img_count += 1
                    last_tsec = tsec
             
            # Time delay.
            while True:
                cpumsec2 = 1000.0*time.perf_counter()
                if (cpumsec2-cpumsec1) >= self.waitmsec:
                    break
            
            # Get next line.
            line = self.lfnt.readline()
        
        # Quit or Done; flush plot events and close trajectory data file.
        
        if self.backend[0:2] == 'WX':
            self.canvas.update()
        elif self.backend[0:2] == 'TK':
            self.canvas.blit(self.ax.bbox)
        elif self.backend[0:2] == 'NB':
            self.canvas.blit(self.ax.bbox)
        elif self.backend == 'MODULE://IPYMPL.BACKEND_NBAGG':
            self.canvas.blit(self.ax.bbox)
        else:  # QT or GTK
            self.canvas.update()
        self.canvas.flush_events()
        
        if self.imgSave and ktot < 0:
            # Save duplicate of last image to ensure final frame in an
            # animated GIF or MP4 video file shows time of intercept.
            sbuff = StringIO()
            sbuff.write("./Ximg/img_%04hd.png" % (img_count))
            img_fname = sbuff.getvalue()
            sbuff.close()
            self.fig.savefig(img_fname, format='png')
                    
        self.lfnt.close()
        self.lfnt = None
