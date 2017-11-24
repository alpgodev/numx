c=======================================================================
c
c     subroutine EXARVO                                    
c
c     Computes Ex-ante Volatility
c
c-----------------------------------------------------------------------
      SUBROUTINE EXARVO (n, cov, w, vol, info)
c      
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       cov    : covariance matrix (n*n)                          double
c       w      : weights (n)                                      double
c
c     OUTPUT 
c       vol    : volatility                                       double
c       info   : diagnostic argument                             integer
c
c     CALL 
c       EXARVA : ex-ante variance
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION vol, w(*), cov(*)
c
c     local variables
      DOUBLE PRECISION var
c
c     external subroutines
      EXTERNAL EXARVA
c     
c     intrinsic functions
      INTRINSIC sqrt
c
c-----------------------------------------------------------------------
c
c     initialization
      vol  = 0.0
      info = 0
c
c     ex-ante variance
      CALL EXARVA ( n, cov, w, var, info)
      IF (info .LT. 0) RETURN
c
c     volatility
      vol = SQRT(var)
c
      RETURN
      END

