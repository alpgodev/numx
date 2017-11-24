c=======================================================================
c
c     subroutine EXPSRA                                     
c
c     Ex-post Sharpe ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPSRA ( n, x, rfr, sharpe, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c          n      : size of x (n > 1)                            integer
c          x      : returns (n)                                   double
c          rfr    : risk-free rate                                double
c
c     OUTPUT 
c          sharpe : Sharpe ratio                                  double
c          info   : diagnostic argument                          integer
c
c     CALL   
c          MV     : computing the mean of a vector
c          VOLAT  : computing the volatility (standard deviation)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION rfr, sharpe
      DOUBLE PRECISION x(*)
c
c     local variables
      DOUBLE PRECISION vol, mret, myzero
      PARAMETER ( myzero = 1.e-30 )
c
c     external subroutines
      EXTERNAL MV, VOLAT  
c
c-----------------------------------------------------------------------
c
c     computing mean return
      CALL MV ( n, x, mret )
c
c     computing volatility
      CALL VOLAT ( n, x, vol, info )
      IF (info .lt. 0) RETURN
c
c     test if volatility is too small
      IF ( vol .le. myzero ) THEN
         info = -104
         RETURN
      ENDIF
c
c     computing sharpe ratio
      sharpe = ( mret - rfr ) / vol
c
      RETURN
      END
c
c=======================================================================
