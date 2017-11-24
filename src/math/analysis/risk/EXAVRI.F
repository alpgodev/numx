c=======================================================================
c
c     subroutine EXAVRI                                      
c
c     Ex-ante Value-at-Risk (upper bound)
c
c-----------------------------------------------------------------------
      SUBROUTINE EXAVRI ( n, cov, rho, w, conflp, valris, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weights (n)                                      double
c       cov    : covariance matrix (n*n)                          double
c       rho    : asset(s) mean return(s) (n)                      double
c       conflp : confidence level parameter ( 0<conflp<1 )        double
c
c     OUTPUT 
c       valris : Value-at-Risk (upper bound)                      double
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
c     argument i/o
      INTEGER n, info
      DOUBLE PRECISION conflp, valris, w(*), cov(*), rho(*)
c
c     local variables
      DOUBLE PRECISION mret, vol, EPS
      PARAMETER ( EPS = 1.E-8 )
c
c     external subroutines
      EXTERNAL EXARET, EXARVO
c     
c     intrinsic functions
      INTRINSIC SQRT, MIN, MAX    
c
c-----------------------------------------------------------------------
c
c     initialization
      valris = 0.0 
      info   = 0
c
c     test confidence probability level (conflp) 
      IF ( (conflp .LT. EPS) .or. (conflp .GT. (1. - EPS)) ) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-ante mean return : rmean
      CALL EXARET ( n, rho, w, mret, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante volatility : vol
      CALL EXARVO ( n, cov, w, vol, info)
      IF (info .LT. 0) RETURN
c
c     computing Value-at-Risk upper bound  
      valris = - vol * SQRT( conflp/(1.-conflp)) + mret
c      valris = MIN(valris, 1.)
c      valris = MAX(valris, -1.)
c
      RETURN
      END
