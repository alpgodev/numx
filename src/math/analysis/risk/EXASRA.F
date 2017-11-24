c=======================================================================
c
c     subroutine EXASRA                                      
c
c     Computes ex-ante Sharpe ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXASRA ( n, cov, rho, w, rfr, sharpe, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weights  (n)                                     double
c       cov    : covariance matrix (n*n)                          double
c       rho    : mean returns (n)                                 double
c       rfr    : risk-free rate                                   double
c
c     OUTPUT 
c       sharpe : Sharpe ratio                                     double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       EXARET  : ex-ante return
c       EXARVO  : ex-ante volatility
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION rfr, sharpe, w(*), cov(*), rho(*)
c
c     local variables
      DOUBLE PRECISION mret, vol, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL EXARET, EXARVO
c
c-----------------------------------------------------------------------
c
c     initialization
      sharpe = 0.0
      info   = 0
c
c     ex-ante mean return
      CALL EXARET ( n, rho, w, mret, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante volatility EXARVO (n, cov, w, vol, info)
      CALL EXARVO ( n, cov, w, vol,  info)
      IF (info .LT. 0) RETURN
c
c     test if volatility is too small
      IF ( vol .LT. EPS ) THEN
         info = -104
         RETURN
      ENDIF
c
c     computing sharpe ratio
      sharpe = ( mret - rfr ) / vol
c
      RETURN
      END
