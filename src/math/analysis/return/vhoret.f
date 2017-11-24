c=======================================================================
c
c     subroutine VHORET                                      
c
c     Arithmetic h-returns for one assets
c
c      r(t+h) = [ price(t+h) - price(t)]/price(t+h) for t = 1,...,n-h
c
c-----------------------------------------------------------------------
      SUBROUTINE vhoret ( n, X, h, R, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n    : number of value(s) (n > h)                  integer
c            X    : asset values (>0) (n)                        double
c            h    : horizon (h > 0)                             integer
c
c     OUTPUT 
c            R    : arithmetic returns  (n-h)                    double
c            info : diagnostic argument                         integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, h, info
      DOUBLE PRECISION X(*), R(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION value, myzero
      PARAMETER ( myzero = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     test if the number of values > h)
      IF (n .le. h) THEN
        info = -2
        RETURN
      ENDIF      
c
c     arithmetic return(s)
      DO i = 1,n-h
         value = X(i)
         IF (value .lt. myzero) THEN
            info = -103
            RETURN
         ENDIF
         R(i) = ( X(i + h)- value ) / value
      ENDDO
c
      RETURN
      END

