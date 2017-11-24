c=======================================================================
c
c     subroutine EXPMVA                                      
c
c     Ex-post Cornish-Fisher (or modified) Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPMVA( n, x, beta, expmv, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c        n       : number of return(s) (n>1)                     integer
c        x       : returns (n)                                    double
c        beta    : quantile of the normal distribution            double
c       
c     OUTPUT :
c        expmv   : ex-post modified Value-at-Risk                 double
c        info    : diagnostic argument                           integer
c
c     CALL :  
c        EXRRET  : ex-post mean return
c        EXPVAR  : ex-post variance    
c        EXPSKE  : ex-post skewness
c        EXPKUR  : ex-post kurtosis
c        dinvnr  : cumulative normal distribution inverse
c        MODVAR  : modified VaR
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o  
      INTEGER n, info    
      DOUBLE PRECISION beta, expmv
      DOUBLE PRECISION x(*)
c
c     local variables      
      DOUBLE PRECISION zc, mexpr, epvar, expsk, expku, p, q, eps
      PARAMETER (eps = 1.E-8)
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL EXPRET, EXPVAR, EXPSKE, EXPKUR, MODVAR
c     
c     intrinsic functions
      INTRINSIC min, max       
c   
c-----------------------------------------------------------------------      
c
c     initialization
      expmv = 0.0
      info  = 0
c
c     test confidence level: beta
      IF ((beta .LT. eps).or.(beta .GT. (1. - eps))) THEN
         info = -102
         RETURN
      ENDIF      
c
c     ex-post mean return
      CALL EXPRET ( n, x, mexpr, info )
      IF (info .lt. 0) RETURN
c
c     ex-post variance
      CALL EXPVAR ( n, x, epvar, info)
      IF (info .lt. 0) RETURN
c
c     ex-post skewness
      CALL EXPSKE ( n, x, expsk, info )
      IF (info .lt. 0) RETURN
c
c     ex-post kurtosis
      CALL EXPKUR ( n, x, expku, info )
      IF (info .lt. 0) RETURN
c
c     computing inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p = beta
      q = 1. - p
      zc = dinvnr(p,q)
c
c     modified VaR (Gaussian)
      CALL MODVAR(mexpr, epvar, expsk, expku, zc, expmv, info )
      IF (info .lt. 0) RETURN      
      expmv = min(expmv, 1.)
      expmv = max(expmv, -1.)
c
      RETURN
      END
