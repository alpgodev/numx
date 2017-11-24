c=======================================================================
c
c     subroutine EXPSOR                                      
c
c     Ex-post Sortino ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPSOR ( n, x, rfr, sortin, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c          n       : size of x (n > 1)                           integer
c          x       : returns (n)                                  double
c          rfr     : risk-free rate                               double
c
c     OUTPUT 
c          sortin  : Sortino ratio                                double
c          info    : diagnostic argument                         integer
c
c     CALL   
c          MV      : mean of a vector
c          SEMVOL  : Computing the semi-volatility
c                  ( volatility(v), v(i) < mean(v(i)) )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/O
      INTEGER n, info
      DOUBLE PRECISION rfr, sortin
      DOUBLE PRECISION x(*)
c
c     local variables  
      DOUBLE PRECISION svol, mean, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL MV, SEMVOL   
c
c-----------------------------------------------------------------------
c
c     initialization
      sortin = 0.0
      info   = 0
c
c     computing mean return
      CALL MV ( n, x, mean )
c
c     computing volatility
      CALL SEMVOL ( n, x, svol, info )
      IF (info .LT. 0) RETURN
c
c     test if volatility is too small
      IF ( svol .LT. EPS ) THEN
         info = -104
         RETURN
      ENDIF
c
c     computing sharpe ratio
      sortin = ( mean - rfr ) / svol
c
      RETURN
      END

