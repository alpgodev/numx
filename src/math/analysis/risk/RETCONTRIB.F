c=======================================================================
c
c     subroutine RETCONTRIB                                  
c
c     Performance (return) Contribution (marginal perf. contribution)
c
c     perf. contribution(i) = w(i)*rho(i) / (Total Return)
c
c-----------------------------------------------------------------------
      SUBROUTINE RETCONTRIB( n, w, rho, mret, mper, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weight(s) (n>=1)                                 double 
c       rho    : asset(s) (mean) return (n)                       double
c
c     OUTPUT 
c       mret   : marginal performance contribution (n)            double
c       mper   : marginal performance contribution in % (n)       double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       PVX    : Vector * Scalar = Vector
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION w(*), rho(*), mret(*), mper(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION EPS, exantr, sum
      PARAMETER (EPS = 1.E-30)
c
c     external subroutines
      EXTERNAL XV, PVX  
c
c-----------------------------------------------------------------------
c
c     initialisation 
      info = 0
      CALL IVX ( n, mret, 0.0 )
      CALL IVX ( n, mper, 0.0 )
c
c--------------------------------------------------------------------      
c
c     ex-ante portfolio mean return := w'*rho
      CALL XV ( n, w, rho, exantr )
c
c     marginal contribution w(i)*rho(i)
      DO i = 1,n
        mret(i) = w(i)*rho(i)
      ENDDO
      sum = 0.
      DO i = 1,n
        sum = sum + mret(i)
      ENDDO
c      
      IF (ABS(sum) .LT. EPS) THEN
        DO i = 1,n
            mper(i) = w(i)
        ENDDO
        info = -1
        RETURN
      ENDIF
      sum = 1./sum
c
c     marginal contribution in % (marginal contribution)  
c     w(i)*rho(i)/(portfolio return)
      CALL PVX ( n, mret, sum, mper )
c
      RETURN
      END
 
