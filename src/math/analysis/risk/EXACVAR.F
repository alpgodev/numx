c=======================================================================
c
c     subroutine EXACVAR                                     
c
c     Ex-ante conditional Gaussian Value-at-Risk (CVaR)
c
c-----------------------------------------------------------------------
      SUBROUTINE EXACVAR ( n, gamma, rho, w, alpha, cvar, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of asset(s)                         integer
c            w      : weight(s) (n)                               double
c            gamma  : covariance matrix (n*n)                     double
c            epskat : Kato sensibility parameter                  double
c            rho    : asset(s) mean return(s) (n)                 double
c            alpha  : probability confidence level                double
c
c     OUTPUT 
c            cvar   : conditional Gaussian Value-at-Risk          double
c            info   : error code                                 integer                      
c
c     CALL 
c           EXARET  : ex-ante return
c           EXARVO  : ex-ante volatility
c           dinvnr  : cumulative normal distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha, cvar
      DOUBLE PRECISION gamma(*), w(*), rho(*)
c
c     local variables
      DOUBLE PRECISION coeff, zc, p, q, meanx, volx, EPS, PI
      PARAMETER ( EPS = 1.E-8, PI = 3.1415926535898 )
c
c     external functions
      DOUBLE PRECISION dinvnr, spmpar
      EXTERNAL dinvnr, spmpar
c
c     external subroutines
      EXTERNAL EXARET, EXARVO
c     
c     intrinsic functions
      INTRINSIC sqrt, exp, min, max      
c
c-----------------------------------------------------------------------
c
c     initialization
      cvar   = 0.0
      info   = 0
c
c     machine precision
c      myzero = spmpar(2)*(1.E+4)
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
      CALL EXARVO ( n, gamma, w, volx, info)
c
c     inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p  = alpha
      q  = 1. - p
      zc = dinvnr(p, q)
c
      coeff = SQRT(2.0*PI)*EXP(((zc*zc)/4.0))*(1.0 - alpha)
      coeff = 1.0/coeff    
c     
c     computing conditional Gaussian Value-at-Risk 
c     CVaR := max[0, min(1, mean - coeff*volatility)] 
      cvar = meanx - volx * coeff
c      cvar = MIN(cvar, 1.)
c      cvar = MAX(cvar, -1.)
c
      RETURN
      END

