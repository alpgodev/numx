c=======================================================================
c
c     subroutine EXPNSHORTFALL                               
c
c     Ex-post normal Shortfall probability
c
c     shortfall := N[(mu - perf.)/vol.]
c 
c                 N[.]  : normal distribution
c                 mu    : threshold return
c                 perf. : portfolio performance
c                 vol.  : portfolio volatility
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            n     : number of return(s) (n > 1)                 integer
c            x     : return(s) (n)                                double
c            mu    : threshold return                             double
c
c     OUTPUT :
c           sfall  : Gaussian shortfall probability               double
c           info   : diagnostic argument                         integer
c
c     CALL   :
c          EXPRET  : ex-post mean return
c          VOLAT   : ex-post volatility
c
c-----------------------------------------------------------------------
c
      SUBROUTINE EXPNSHORTFALL ( n, x, mu, sfall, info)
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION mu, sfall
      DOUBLE PRECISION x(*)
c
c     local variables      
      DOUBLE PRECISION u, volx, meanx, pi
      PARAMETER ( pi = 3.1415926535898 )
c
c     external subroutines
      EXTERNAL EXPRET, VOLAT
c     
c     intrinsic functions
      INTRINSIC sqrt, exp, min, max      
c
c-----------------------------------------------------------------------
c
c     initialization
      sfall = 0.0
      info  = 0
c
c     ex-post mean return
      CALL EXPRET ( n, x, meanx, info )
      IF (info .lt. 0) RETURN
c
c     ex-post volatility
      CALL VOLAT ( n, x, volx, info )
      IF (info .lt. 0) RETURN
c     
c     computing Normal Shortfall probability 
      u     = (mu - meanx)/volx
      sfall = exp(-(u*u)/2)/sqrt(2.0*pi) 
      sfall = min(sfall, 1.)
      sfall = max(sfall, -1.)
c
      RETURN
      END
