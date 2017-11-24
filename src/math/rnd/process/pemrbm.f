c=======================================================================
c
c     subroutine PEMRBM                                       
c
c     Mean Reversion process (Arithmetic Ornstein-Uhlenbeck)
c
c           dS(t) = -theta*(S(t)-mu)dt + sigma*dW(t)
c
c     Method: Euler discretization of SDE
c
c                S(1)   = S
c                a      = exp(- delta*dt )
c                b      = sqrt((1 - exp(- 2*delta*dt))/(2*delta))
c                S(t+1) = S(t)*a + mu*(1 - a) + sigma*b*N(0,1)
c
c                with dt = T/(N-1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE pemrbm ( T, N, S, theta, mu, sigma, y, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            S      : initial value                               double
c            theta  : velocity of the reversion                   double
c            mu     : mean reversion                              double
c            sigma  : process volatility                          double
c
c     OUTPUT :
c            y      : generated step (N)                          double
c            info   : diagnostic information                     integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, T, info
      DOUBLE PRECISION S, theta, mu, sigma, y(*)
c
c     local variables
      INTEGER i, M
      DOUBLE PRECISION dt, a, b
c
c     external functions
      DOUBLE PRECISION rn
      EXTERNAL rn
c
c     intrinsic function
      INTRINSIC DFLOAT, EXP, SQRT
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
c     coefficients
      a = EXP( - theta * dt)
      b = SQRT((1. - EXP(- 2.* theta * dt))/(2. * theta))
c
c     mean-reversion process 
      y(1) = S
      DO i = 1,M
         y(i + 1) = (a*y(i)) + (mu*(1. - a)) + (b*sigma*rn()) 
      ENDDO
c
      RETURN
      END

