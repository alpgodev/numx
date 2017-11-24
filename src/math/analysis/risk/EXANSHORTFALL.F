c=======================================================================
c
c     subroutine EXANSHORTFALL                              
c
c     Ex-ante Gaussian Shortfall probability
c
c-----------------------------------------------------------------------
      SUBROUTINE EXANSHORTFALL ( n, cov, rho, w, mu, sfall, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weight(s) (n)                                    double
c       cov    : covariance matrix (n*n)                          double
c       rho    : asset(s) mean return(s) (n)                      double
c       mu     : threshold return                                 double
c
c     OUTPUT 
c       sfall  : Gaussian shortfall probability                   double
c       info   : error code                                      integer  
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
      DOUBLE PRECISION mu, sfall, cov(*), w(*), rho(*)
c
c     local variables
      DOUBLE PRECISION u, meanx, volx, PI
      PARAMETER ( PI = 3.1415926535898 )
c
c     external subroutines
      EXTERNAL EXARET, EXARVO
c     
c     intrinsic functions
      INTRINSIC exp, sqrt, min, max    
c
c-----------------------------------------------------------------------
c
c     initialization
      sfall = 0.0
      info  = 0
c
c     ex-ante mean 
      CALL EXARET ( n, rho, w, meanx, info)
c      
c     ex-ante volatility
      CALL EXARVO ( n, cov, w, volx, info)
c     
c     Normal Shortfall probability 
      u = (mu - meanx)/volx
      sfall = EXP(- 0.5*(u*u))/SQRT(2.0*PI) 
      sfall = MIN(sfall, 1.)
      sfall = MAX(sfall, -1.)
c
      RETURN
      END
