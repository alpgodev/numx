c=======================================================================
c
c     subroutine PEGARCH                                     
c
c     Autoregressive Conditionally Heteroskedastic (ARCH) process generator
c
c                e(t) = sigma(t)*u(t)
c
c                where
c
c                sigma(t) = sqrt[w + alpha(1)*e(t-1)^2 + ... + alpha(p)*e(t-p)^2]
c
c-----------------------------------------------------------------------
      SUBROUTINE pegarch ( T, N, P, Q, w, alpha, beta, y, sigma, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            T      : time (days, months, years)                 integer
c            N      : number of sub-division(s)                  integer 
c            P      : AR process order (P >= 0)                  integer
c            Q      : MA process order (Q >= 0)                  integer
c            w      : constant parameter                          double
c            alpha  : GARCH P-parameter, vector (P)               double
c            beta   : GARCH Q-parameter, vector (Q)               double
c
c     OUTPUT 
c            y      : generated process (N)                       double
c            sigma  : generated conditional volatility (N)        double
c            info   : diagnostic information                     integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER N, T, P, Q, info
      DOUBLE PRECISION w, alpha(*), beta(*), y(*), sigma(*)
c
c     local variables 
      INTEGER M, i, j, k
      DOUBLE PRECISION sum, epsilon
      PARAMETER (epsilon = 1.E-15)     
c
c     external subroutines
      EXTERNAL PEARCH  
c
c     external functions
      DOUBLE PRECISION rn
      EXTERNAL rn
c     
c     intrinsic functions
      INTRINSIC max, sqrt  
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     maximum order max(P,Q)
      M = max(P,Q)
c
c     case Q = 0 (ARCH process)
      IF (Q .EQ. 0) THEN
        CALL pearch ( T, N, P, w, alpha, y, sigma, info )
        RETURN
      ENDIF      
c
c     test if M < 0 or M > N
      IF ((M .LT. 1).OR.(M .GT. N)) THEN
        info = -3101
        RETURN
      ENDIF
c
c     initialization step 
      sum = 0.
      DO j = 1,P
        sum = sum + alpha(j)
      ENDDO
      DO j = 1,Q
        sum = sum + beta(j)
      ENDDO
c
c     test if sum(alpha+beta) < 1. 
      IF (sum .GE. (1. - epsilon)) THEN
        info = -3102
        RETURN
      ENDIF
c      
c     loop
      DO i = 1,M
        sigma(i) = sqrt(w/(1.-sum))
        y(i)     = sigma(i)*rn()
      ENDDO  
c
c     loop      
      DO i = 1,N-M
c
c        conditional volatility    
         sum = w   
         IF (P .GE. 1) THEN 
            DO j = 1,P
              k = i + P - j
              sum = sum + alpha(j)*y(k)*y(k)
            ENDDO
         ENDIF   
         IF (Q .GE. 1) THEN
            DO j = 1,Q
                k = i + Q - j
                sum = sum + beta(j)*sigma(k)*sigma(k)
            ENDDO
         ENDIF   
         sigma(i + M) = sqrt(sum) 
c
c        GARCH(P,Q) process
         y(i + M) = sigma(i + M)*rn()
      ENDDO
c
      RETURN
      END
