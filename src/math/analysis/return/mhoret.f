c=======================================================================
c
c     subroutine MHORET                                     
c
c     Arithmetic returns for prices matrix
c
c       r(t+h) = [ price(t+h) - price(t)]/price(t+h) for t = 1,...,n-h
c
c----------------------------------------------------------------------
      SUBROUTINE mhoret ( n, p, X, h, R, info )
c----------------------------------------------------------------------
c
c     INPUT :
c            n     : number of prices (n > h)                   integer
c            p     : number of asset(s) (p >= 1)                integer
c            X     : assets value (>0) (n*p)                     double
c            h     : horizon (h > 0)                            integer
c
c     OUTPUT :
c            R     : aritmetic returns ( n-h * p)                double    
c            info  : diagnostic argument                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, h, info
      DOUBLE PRECISION X(n,*), R(n-h,*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION value, myzero
      PARAMETER ( myzero = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      value = 0.0
c
c     test if the number of values > h
      IF (n .le. h) THEN
        info = -2
        RETURN
      ENDIF
c
c     artithmetic return(s)
      DO j = 1,p
        DO i = 1,n - h
           value = X(i, j)
           IF (value .lt. myzero) THEN
              info = -103
              RETURN
           ENDIF
           R(i, j) = ( X(i + h, j) - value )/value
         ENDDO
      ENDDO
c
      RETURN
      END

