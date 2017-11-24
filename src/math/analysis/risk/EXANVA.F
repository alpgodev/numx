c=======================================================================
c
c     subroutine EXANVA                                      
c
c     Ex-ante Gaussian Value-at-Risk (parametric method)
c
c     VaR(alpha) := w'*rho - z(alpha)*SQRT(w'*Cov*w)
c
c-----------------------------------------------------------------------
      SUBROUTINE EXANVA ( n, cov, rho, w, alpha, nvar, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weight(s) (n)                                    double
c       cov    : covariance matrix (n*n)                          double
c       rho    : asset(s) mean return(s) (n)                      double
c       alpha  : probability confidence level                     double
c
c     OUTPUT 
c       nvar   : Gaussian Value-at-Risk                           double
c       info   : error code                                      integer                 
c
c     CALL   
c       EXARET  : ex-ante return
c       EXARVO  : ex-ante volatility
c       dinvnr  : cumulative normal distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha, nvar, cov(*), w(*), rho(*)
c
c     local variables
      DOUBLE PRECISION zc, p, q, meanx, volx, EPS
      PARAMETER ( EPS = 1.e-8 )
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL EXARET, EXARVO
c     
c     intrinsic functions
      INTRINSIC min, max   
c
c-----------------------------------------------------------------------
c
c     initialization
      nvar = 0.0
      info = 0
c
c     test confidence level
      IF ( (alpha .LT. EPS) .or. (alpha .GT. (1. - EPS)) ) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-ante mean 
      CALL EXARET ( n, rho, w, meanx, info)
c      
c     ex-ante volatility
      CALL EXARVO ( n, cov, w, volx, info)
c
c     inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p = alpha
      q = 1. - p
      zc = dinvnr(p,q)
c     
c     computing Normal VaR := max[0, min(1, -zc*volatility + mean)] 
      nvar = - volx * zc + meanx
c      nvar = MIN(nvar, 1.)
c      nvar = MAX(nvar, -1.)
c
      RETURN
      END
