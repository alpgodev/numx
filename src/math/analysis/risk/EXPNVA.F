c=======================================================================
c
c     subroutine EXPNVA                                      
c
c     Ex-post normal (Gaussian) Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPNVA ( n, x, alpha, nvar, info)
c-----------------------------------------------------------------------
c
c     INPUT :
c            n     : number of returns                           integer
c            x     : returns (n)                                  double
c            alpha : probability level (0 < alpha < 1)            double
c
c     OUTPUT :
c            nvar  : Gaussian Value-at-Risk                       double
c            info  : = 0 successful exit                         integer
c
c     CALL   :
c          EXPRET  : ex-post mean return
c          VOLAT   : ex-post volatility
c          dinvnr  : cumulative normal distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha, nvar
      DOUBLE PRECISION x(*)
c
c     local variables
      DOUBLE PRECISION zc, p, q, volx, meanx, eps
      PARAMETER ( eps = 1.E-8 )
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL EXPRET, VOLAT
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
c     test confidence level: alpha
      IF ((alpha .LT. eps).or.(alpha .GT. (1.-eps))) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-post mean return
      CALL EXPRET ( n, x, meanx, info )
      IF (info .lt. 0) RETURN
c
c     ex-post volatility
      CALL VOLAT ( n, x, volx, info )
      IF (info .lt. 0) RETURN
c
c     computing inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p = alpha
      q = 1. - p
      zc = dinvnr(p,q)
c     
c     computing Normal -VaR := max[0, min(1, -zc*volatility - mean)] 
      nvar = - volx * zc + meanx 
      nvar = min(nvar, 1.)
      nvar = max(nvar, -1.)
c
      RETURN
      END
