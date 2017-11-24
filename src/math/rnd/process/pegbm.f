c=======================================================================
c
c     subroutine PEGBM                                      
c
c     Simulation of a Geometric Brownian Motion
c
c               dSt = St*[mu*dt + sigma*dWt]
c
c     Method: Euler discretization of SDE
c
c                S(1)    = S
c                S(1+t)  = S(t)*[1 + mu*dt +  sigma*sqrt(dt)*N(0,1)]
c
c                with dt = T/(N-1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE pegbm ( T, N, S, mu, sigma, y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            S      : initial value                               double
c            mu     : process drifts                              double
c            sigma  : process volatility                          double
c
c     OUTPUT 
c            y      : generated step (N)                          double
c            info   : diagnostic information                     integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, T, info
      DOUBLE PRECISION S, mu, sigma, y(*)
c
c     local variables
      INTEGER i, M
      DOUBLE PRECISION dt
c
c     external functions
      DOUBLE PRECISION rn
      EXTERNAL rn
c
c     intrinsic function
      INTRINSIC DFLOAT, SQRT
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     discetization step 
      M = N - 1
      dt = (DFLOAT(T)/M)
c
c     standard Brownian motion
      y(1) = S
      DO i = 1,M
         y(i + 1) = y(i)*(1.0 + mu*dt + sigma*SQRT(dt)*rn()) 
      ENDDO
c
      RETURN
      END
