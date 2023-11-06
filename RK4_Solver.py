# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name

# File: RK4_Solver.py
# Auth: G. E. Deschaines
# Date: 19 May 2015
# Prog: ODE solver using Runge-Kutta 4th order integration method.
# Desc: Applies the Runge-Kutta 4th order (RK4) integration method to
#       solve a system of first Order Differential Equations (ODEs) of
#       the form dS[i] = dotS(i,S[i]), such that each element of the
#       state vector S are calculated as the weighted sum of four
#       approximations of dS[i] and S[i] = S0[i] + dS[i]*h for each
#       i from 0 to n-1.  The S[0] element of the state vector holds
#       incremented time and the associated state derivative dS[0]
#       value must be specified as 1.
#
# Disclaimer:
#
# See DISCLAIMER

import sys

try:
    import numpy as np
except ImportError:
    print("* Error: NumPy package required.")
    print("         Suggest installing the SciPy stack.")
    sys.exit()
  
class RK4_Solver:
    """
    Runge-Kutta 4th-Order Solver
    """
    ## Constructor
  
    def __init__(self, h, n):
        """
        Instantation initializer, where:
            h = integration step size
            n = number of state variables
        """
        self.S      = np.zeros(n)  # temp storage for states
        self.Sinit  = np.zeros(n)  # save of initial states
        self.Sprev  = np.zeros(n)  # save of previous states
        self.Scurr  = np.zeros(n)  # save of current states
        self.dS     = np.zeros(n)  # temp storage for state derivatives 
        self.n      = n            # number of state variables
        self.h      = h            # integration step size 
        self.hh     = 0.5*h        # integration half step size
        self.sixth  = 1.0/6.0      # estimated derivatives weighting factor
        self.h6th   = h*self.sixth 

    
    ## Private Methods

    def _substep(self, h, St, dS):
        """
        Returns vector S which holds the solution to the equation
        S[i] = St[i] + dS[i]*h for i from 0 to n-1, where:
            h  = integration step size
            n  = number of state variables
            St = state vector at t (i.e., [t,x,v])
            dS = state 1st derivatives vector at t (i.e., [1,dx/dt,dv/dt]).
        """
        S = np.zeros(self.n)

        for i in range(self.n) :
            S[i] = St[i] + dS[i]*h

        return S
        
  
    ## Public Methods
  
    def init(self, S) :
        """
        Initializes the state vectors.
        """
        self.Sinit = S
        self.Sprev = S
        self.Scurr = S
        
        
    def step(self, Scurr, dotS):
        """
        Returns vector S which holds the RK4 solution to the
        equation S = Scurr + dotS*h, where:
            Scurr = state vector at current time (i.e., [t,x,v])
            dotS  = function of the form dotS(n,S) containing
                    system of 1st order differential equations
                    to integrate over the time step h.
        """ 
        K1 = dotS(self.n, Scurr)
        K2 = dotS(self.n, self._substep(self.hh, Scurr, K1))
        K3 = dotS(self.n, self._substep(self.hh, Scurr, K2))
        K4 = dotS(self.n, self._substep(self.h,  Scurr, K3))

        for i in range(self.n) :
            self.S[i] = Scurr[i] + \
                        (K1[i] + 2.0*(K2[i] + K3[i]) + K4[i])*self.h6th

        self.Sprev = self.Scurr.copy()
        self.Scurr = self.S.copy()
        
        return self.S
        
    
    def get_dSprev(self, dotS):
        """
        Returns vector dS which holds the RK4 weighted average
        values of the state derivatives computed from the previous 
        state and the given state derivatives function dotS.

        """
        K1 = dotS(self.n, self.Sprev)
        K2 = dotS(self.n, self._substep(self.hh, self.Sprev, K1))
        K3 = dotS(self.n, self._substep(self.hh, self.Sprev, K2))
        K4 = dotS(self.n, self._substep(self.h,  self.Sprev, K3))

        for i in range(self.n) :
            self.dS[i]  = (K1[i] + 2.0*(K2[i] + K3[i]) + K4[i])*self.sixth

        return self.dS
        
        
    def get_dSnext(self, dotS):
        """
        Returns vector dS which holds the RK4 weighted average
        values of the state derivatives computed from the current 
        state and the given state derivatives function dotS.

        """
        K1 = dotS(self.n, self.Scurr)
        K2 = dotS(self.n, self._substep(self.hh, self.Scurr, K1))
        K3 = dotS(self.n, self._substep(self.hh, self.Scurr, K2))
        K4 = dotS(self.n, self._substep(self.h,  self.Scurr, K3))

        for i in range(self.n) :
            self.dS[i] = (K1[i] + 2.0*(K2[i] + K3[i]) + K4[i])*self.sixth

        return self.dS
        
        
    def get_Sinit(self) :
        """
        Returns initial state vector
        """
        return self.Sinit
        
    
    def get_Sprev(self) :
        """
        Returns state vector for the previous time step.
        """
        return self.Sprev
        
    
    def get_Scurr(self) :
        """
        Returns state vector for the current time step.
        """
        return self.Scurr
    