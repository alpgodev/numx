c=======================================================================
c
c     subroutine PEBM
c
c     Simulation of a Brownian process (Wiener process)
c 
c                dW(t) = W(t+dt) - W(t) := N(0,1)*std
c
c     Method: Euler discretization of SDE
c
c                w(1)    = S
c                W(1+t)  = W(t) + sqrt(dt)*N(0,1) for t in [1,T]
c
c                with dt = T/(N-1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE pebm ( T, N, S, y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            S      : initial value                               double
c
c     OUTPUT 
c            y      : generated process (N)                       double
c            info   : diagnostic information                     integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER T, N, info
      DOUBLE PRECISION S, y(*)
c
c     local variables
      INTEGER i, M
      DOUBLE PRECISION dt
c
c     external functions
      DOUBLE PRECISION rn           
      EXTERNAL rn
c      
c     intrinsic functions 
      INTRINSIC dfloat, sqrt   
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     discetization step  
      M = N - 1
      dt = (dfloat(T)/M)
c
c     Wienner process - Brownian motion
      y(1) = S
      DO i = 1,M
         y(i + 1) = y(i) + sqrt(dt) * rn() 
      ENDDO
c
      RETURN
      END
