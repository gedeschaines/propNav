# -*- coding: utf-8 -*-

# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name
# pylint: disable=anomalous-backslash-in-string,bad-continuation
# pylint: disable=multiple-statements,redefined-outer-name,global-statement

"""
FILE:  clipLib.py
DATE:  06 DEC 2023
AUTH:  G. E. Deschaines
DESC:  Methods to determine clipping of a given line segment with
       the edges of a viewing pyramidal frustum.  These methods
       were derived from algorithms presented on pages 152-155 in
       Chapter 3 of "Procedural Elements for Computer Graphics"
       by David F. Rogers, published by McGraw-Hill, Inc., 1985.
REFS:
     
  [1] This Python script was refactored from cliplib.c available at:
      https://github.com/gedeschaines/threeD/blob/master/src/cliblib.c


Disclaimer:  See DISCLAIMER file.

"""

from draw3D import DBG_LVL, Pnt3D

mxvcnt = 32       # maximum vertices in clipped polygon 
zmin   = 0.1      # minimum z clipping distance
zmax   = 10000.0  # maximum z clipping distance


def edgeCode(edge, a_pt):
    """
    Calculates edge code for given pyramidal frustum edge and
    polygon vertex point.
    """

    # Load point into viewing pyramid space.

    x =  a_pt.Y
    y = -a_pt.Z
    z =  a_pt.X

    # Initialize edge code.

    code = 0

    # Calculated edge code.
    
    def switch_cases(edge):
        
        def case_1(x, y, z, code):
            if ( x == -z ): return code  # on left edge
            if ( x <  -z ): code = -1    # outside left edge
            else          : code =  1    # inside left edge
            return code
        
        def case_2(x, y, z, code):
            if ( x ==  z ): return code  # on right edge
            if ( x >   z ): code = -2    # outside right edge
            else          : code =  2    # inside right edge
            return code
        
        def case_3(x, y, z, code): 
            if ( y == -z ): return code  # on bottom edge
            if ( y <  -z ): code = -4    # below bottom edge
            else          : code =  4    # above bottom edge 
            return code
        
        def case_4(x, y, z, code):
            if ( y ==  z ): return code  # on top edge
            if ( y >   z ): code = -8    # above top edge
            else          : code =  8    # below top edge
            return code
        
        def case_5(x, y, z, code): 
            if ( z ==  zmax ): return code  # on zmax
            if ( z >   zmax ): code = -16   # in front of zmax
            else             : code =  16   # behind zmax     
            return code
        
        def case_6(x, y, z, code): 
            if ( z ==  zmin ): return code  # on zmin
            if ( z <   zmin ): code = -32   # behind zmin
            else             : code =  32   # in front of zmin
            return code
        
        cases = {1:case_1, 2:case_2, 3:case_3, 4:case_4, 5:case_5, 6:case_6}
       
        return cases[edge]

    switch_case = switch_cases(edge)

    return switch_case(x, y, z, code)


def edgeClip(edge, pt_s, pt_e):
    """
    Determines pyramidal frustum edge clipping of given line segment.
    """
    
    xs =  pt_s.Y
    ys = -pt_s.Z
    zs =  pt_s.X
    xe =  pt_e.Y
    ye = -pt_e.Z
    ze =  pt_e.X
   
    def switch_cases(edge):
       
        def case_1(xs, ys, zs, xe, ye, ze): 
            # left edge intercept
            k   = xe-xs
            t   = (zs+xs)/(zs-ze-k)
            xsp = k*t + xs
            ysp = (ye-ys)*t + ys
            zsp = -xsp
            return xsp, ysp, zsp
       
        def case_2(xs, ys, zs, xe, ye, ze):
            # right edge intercept
            k   = xe-xs
            t   = (zs-xs)/(zs-ze+k)
            xsp = k*t + xs
            ysp = (ye-ys)*t + ys
            zsp = xsp
            return xsp, ysp, zsp
       
        def case_3(xs, ys, zs, xe, ye, ze):
            # bottom edge intercept
            k   = ye-ys
            t   = (zs+ys)/(zs-ze-k)
            xsp = (xe-xs)*t + xs
            ysp = k*t + ys
            zsp = -ysp
            return xsp, ysp, zsp
       
        def case_4(xs, ys, zs, xe, ye, ze):
            # top edge intercept
            k   = ye-ys
            t   = (zs-ys)/(zs-ze+k)
            xsp = (xe-xs)*t + xs
            ysp = k*t + ys
            zsp = ysp
            return xsp, ysp, zsp
       
        def case_5(xs, ys, zs, xe, ye, ze):
            # max z clip plane intercept
            k   = ze-zs
            t   = (zmax-zs)/k
            xsp = (xe-xs)*t + xs
            ysp = (ye-ys)*t + ys
            zsp = zmax
            return xsp, ysp, zsp
       
        def case_6(xs, ys, zs, xe, ye, ze):
            # min z clip plane intercept
            k   = ze-zs
            t   = (zmin-zs)/k
            xsp = (xe-xs)*t + xs
            ysp = (ye-ys)*t + ys
            zsp = zmin
            return xsp, ysp, zsp
       
        cases = {1:case_1, 2:case_2, 3:case_3, 4:case_4, 5:case_5, 6:case_6}
       
        return cases[edge]

    switch_case = switch_cases(edge)
    
    xsp, ysp, zsp = switch_case(xs, ys, zs, xe, ye, ze)
    
    pt_i = Pnt3D(zsp, xsp, -ysp)
    
    return pt_i


def polyClip(pcnt, vcnt, vlist):
    """
    Clips given polygon to 3D viewing pyramidal frustum.
    """

    while True:
        # Check all polygon points against each frustum edge.
        pcntp1 = pcnt + 1
        jcnt   = 0
        pt_S   = vlist[pcnt][1]
        cs     = edgeCode(pcnt, pt_S)
        if ( cs >= 0 ):
            # pt_S inside or on frustum edge - save.
            jcnt                = jcnt + 1
            vlist[pcntp1][jcnt] = pt_S
      
        for icnt in range(2, vcnt[pcnt]+1):
            # Check all subsequent points along polygon.
            pt_E = vlist[pcnt][icnt]
            ce   = edgeCode(pcnt, pt_E)
            if cs != ce:
                # Line segment intercepts frustum edge.
                if cs < ce:
                    # pt_S left of pt_E.
                    pt_X = edgeClip(pcnt, pt_S, pt_E)
                else:
                    # pt_E left of pt_S.
                    pt_X = edgeClip(pcnt, pt_E, pt_S)
                # Save this intercept.
                jcnt                = jcnt + 1
                vlist[pcntp1][jcnt] = pt_X
            if icnt < vcnt[pcnt]:
                # Not last polygon point.
                pt_S = pt_E
                cs   = ce
                if cs >= 0:
                    # Save this point.
                    jcnt                = jcnt + 1
                    vlist[pcntp1][jcnt] = pt_S

        if jcnt > 0:
            # Close polygon.
            jcnt                = jcnt + 1
            vlist[pcntp1][jcnt] = vlist[pcntp1][1]
     
        pcnt       = pcntp1
        vcnt[pcnt] = jcnt
        if DBG_LVL > 5:
            for i in range(1, vcnt[pcnt]+1):
                print("polyClip:  %d %d %f %f %f" % \
                      (pcnt, i, vlist[pcnt][i].X, vlist[pcnt][i].Y, vlist[pcnt][i].Z))
        
        if ( pcnt == 7 ) or ( jcnt == 0 ):
            break

    return pcnt
