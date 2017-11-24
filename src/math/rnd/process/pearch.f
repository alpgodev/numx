c=======================================================================
c
c     subroutine PEARCH                                     
c
c     Autoregressive Conditionally Heteroskedastic (ARCH) process generator
c
c       e(t) = sigma(t)*u(t)
c
c       where
c
c       sigma(t) = sqrt[w + alpha(1)*e(t-1)^2 + ... + alpha(p)*e(t-p)^2]
c
c-----------------------------------------------------------------------
      SUBROUTINE pearch ( T, N, P, w, alpha, y, sigma, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            P      : process order (P > 0)                      integer
c            w      : constant parameter                          double
c            alpha  : ARCH parameter, vector (P)                  double
c
c     OUTPUT :
c            y      : generated process (N)                       double
c            sigma  : generated conditional variance (N)          double
c            info   : diagnostic information                     integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, T, P, info
      DOUBLE PRECISION w, alpha(*), y(*), sigma(*)
c
c     local variables 
      INTEGER i, j, k
      DOUBLE PRECISION sum
c
c     external functions
      DOUBLE PRECISION rn
      EXTERNAL rn
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     test time
      IF (T .LT. 0) RETURN
c
c     test if P < 0 or P > N
      IF ((P .LT. 1).OR.(P .GT. N)) THEN
        info = -3101
        RETURN
      ENDIF
c
c     initialization step
      sum = 0.
      DO j = 1,P
        sum = sum + alpha(j)
      ENDDO
c
c     test if sum(alpha) < 1. 
      IF (sum .GE. 1.) THEN
        info = -3102
        RETURN
      ENDIF
c      
c     initialization of initial value(s)
c     y(1),...,y(p) and volatility 
      DO i = 1,P
        sigma(i) = sqrt(w/(1.-sum))
        y(i)     = sigma(i)*rn()
      ENDDO  
c
c     loop      
      DO i = 1,N-P
c
c        conditional volatility    
         sum = w    
         DO j = 1,P
            k = i + P - j
            sum = sum + alpha(j)*y(k)*y(k)
         ENDDO
         sigma(i + P) = sqrt(sum) 
c
c        ARCH(P) process   
         y(i + P) = sigma(i + P)*rn()
      ENDDO
c
      RETURN
      END
