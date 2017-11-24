c=======================================================================
c
c     subroutine RETLACK
c
c     Arithmetic returns form prices matrix
c
c       r(t+h) = [ price(t+h) - price(t)]/price(t+h) for t = 1,...,n-h
c
c     Take into account of missing data (lack of data), 
c     with historical data structure
c
c     X(n,p) = | x(1,1), ..., x(1,p) |  
c              | x(2,1), ..., -1000  |
c              | x(3,1), ..., x(3,p) |
c              |   .   , ...,   .    |
c              |   .   , ...,   .    |
c              | x(n,1), ..., x(n,p) |
c
c     where x(2,p) = -1000 denote a missing return
c
c----------------------------------------------------------------------
      SUBROUTINE retlack ( n, p, X, h, R, info )
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
      DOUBLE PRECISION missing, price, eps
      PARAMETER ( eps = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      missing = -1000
c
c     test if the number of values > h
      IF (n .LE. h) THEN
        info = -2
        RETURN
      ENDIF
c
c     artithmetic return(s)
      DO j = 1,p
        DO i = 1,n - h
           price = X(i, j)
           IF ((price .LT. eps) .AND. (price .GT. missing)) THEN
              info = -103
              RETURN
           ENDIF
           IF ((X(i + h,j) .GT. missing).AND.(price .GT. missing))
     &     THEN
            R(i, j) = ( X(i + h, j) - price )/price
           ELSE
            R(i, j) = -1000
           ENDIF
         ENDDO
      ENDDO
c
      RETURN
      END

