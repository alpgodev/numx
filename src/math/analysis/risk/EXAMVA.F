c=======================================================================
c
c     subroutine EXAMVA                                      
c
c     Ex-ante Cornish-Fisher (or modified) Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXAMVA(n, p, w, x, cov, alpha, dwork, examv, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n     : number of values                            integer
c            p     : number of assets (p>=1)                     integer
c            w     : weights (p)                                  double
c            x     : asset(s) values (n*p)                        double
c            cov   : covariance matrix (p*p)                      double
c            alpha : probability confidence level                 double
c
c     WORKSPACE 
c            dwork : p*( 6*p + n + 2 )                            double
c
c     OUTPUT 
c            examv : ex-ante Cornish-Fisher VaR                   double
c            info  : diagnostic argument                         integer
c
c     CALL   
c          EXARET  : ex-ante mean return
c          EXARVA  : ex-ante variance
c          EXARSK  : ex-ante skewness
c          EXARKU  : ex-ante kurtosis
c          dinvnr  : cumulative normal distribution inverse
c          MODVAR  : modified Value-at-Risk
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION alpha, examv
      DOUBLE PRECISION w(*), x(*), cov(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables      
      INTEGER pdrho, pdcov, pdwcov, pdskew
      DOUBLE PRECISION exarho, robvar, exaske, exakur, zc, q, EPS
      PARAMETER ( EPS = 1.E-8 )
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL COVMVM, EXARET, EXARVA, EXARSK, EXARKU, MODVAR
c     
c     intrinsic functions
      INTRINSIC MIN, MAX       
c
c-----------------------------------------------------------------------
c
c     initialization 
      examv = 0.0
c
c     pointers for double work space : dwork
c     -------------------------------------- 
      pdrho  = 1
c     pdrho : pointer for average return, so ( p ) more
      pdcov  = pdrho + ( p )
c     pdcov : pointer for cov matrix, so ( p*p ) more     
      pdwcov = pdcov + ( p*p)
c     pdwcov : pointer dwork COVMVM, so ( p*n ) more       
      pdskew = pdwcov + ( p*n )
c     pdskew : pointer dwork EXARSK needs p*( 5*p + 1 )
c
c     Total dwork = p*( 6*p + n + 2 )  
c
c----------------------------------------------------------------------
c
c     test confidence level: alpha
      IF ((alpha .LT. EPS).or.(alpha .GT. (1. - EPS))) THEN
         info = -102
         RETURN
      ENDIF      
c
c     covariance matrix and average return 
      CALL COVMVM ( n, p, x, dwork(pdwcov), dwork(pdrho), 
     &            dwork(pdcov), info)
      IF (info .LT. 0) RETURN
c
c     ex-ante mean return
      CALL EXARET ( p, dwork(pdrho), w, exarho, info )
      IF (info .LT. 0) RETURN
c
c     ex-ante variance
      CALL EXARVA ( p, dwork(pdcov), w, robvar, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante skewness
      CALL EXARSK ( n, p, w, x, cov, dwork(pdskew), exaske, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante kurtosis
      CALL EXARKU ( n, p, w, x, cov, exakur, info)
      IF (info .LT. 0) RETURN
c
c     inverse of the N(0,1) cumulative distribution
c     call dinvnr(p,q)
      q  = 1. - alpha
      zc = dinvnr(alpha, q)
c
c     Cornish-Fisher (modified) Value-at-Risk     
      CALL MODVAR(exarho, robvar, exaske, exakur, zc, examv, info)
      IF (info .LT. 0) RETURN
c 
c     ensure mod. VaR in [-1, 1]
c      examv = MIN(examv, 1.)
c      examv = MAX(examv, -1.)
c
      RETURN
      END
