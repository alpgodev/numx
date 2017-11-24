c=======================================================================
c
c     subroutine PEMRJBM
c
c     Mean Reversion with Jump process (Marlim Model)
c
c                S(t) = delta*(mu - S(t))dt + sigma*dW(t) + dq
c
c                Poisson process dq
c
c                dq = 0   with probability (1-lambda)dt
c                   = phi with probability lambda*dt
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
c
c     INPUT :
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            S      : initial value                               double
c            delta  : velocity of the reversion                   double
c            mu     : mean reversion                              double
c            sigma  : process volatility                          double
c            lambda : jump frequency                               double
c            phi    : jump level                                  double
c
c     OUTPUT :
c            y      : generated step (N)                          double
c            info   : diagnostic information                     integer 
c
c-----------------------------------------------------------------------
c
      SUBROUTINE pemrjbm ( T, N, S, delta, mu, sigma, lambda, phi,
     &                     y, info )
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, T, info
      DOUBLE PRECISION S, delta, mu, sigma, lambda, phi, y(*)
c
c     local variables 
      INTEGER i, M
      DOUBLE PRECISION dt, jump, u, v, a, b
c
c     external functions
      DOUBLE PRECISION rn, ranf
      EXTERNAL rn, ranf
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      jump = 0 
c
c     discetization step 
      M  = N - 1
      IF (M .EQ. 0) THEN
        info = -1
        RETURN
      ENDIF
      dt = (dfloat(T)/M)
c
c     local parameters
      a = exp( - delta * dt)
      b = sqrt((1. - exp(- 2.* delta * dt))/(2. * delta))
      v = lambda * dt
c      
c     mean-reversion with jump
      y(1) = S
      DO i = 1,M
c
c        uniform random variable U[0,1]      
         u = ranf()
         IF (u .lt. v) THEN
            jump = phi
         ELSE
            jump = 0
         ENDIF
         y(i + 1) = (a*y(i)) + (mu*(1. - a)) + (b*sigma*rn()) + jump
      ENDDO
c
      RETURN
      END
