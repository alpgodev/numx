c=======================================================================
c
c     subroutine VARCONTRIB                                  
c
c     Value-at-Risk Contribution (marginal risk contribution)
c
c     vol(i) = -w(i)*z(alpha)*[w'*Cov] / (Total Value-at-Risk)
c
c-----------------------------------------------------------------------
      SUBROUTINE VARCONTRIB( n, cov,  rho, w, alpha, ivar, ipvar, info)
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
c       mvar   : marginal Value-at-Risk contribution (n)          double
c       ivar   : incremental Value-at-Risk contribution (n)       double
c       ipvar  : incremental Value-at-Risk contribution in % (n)  double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       PMV    : Matrix * Vector = Vector
c       PVX    : Vector * Scalar = Vector
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha
      DOUBLE PRECISION w(*), cov(*), rho(*), ivar(*), ipvar(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION vol, zc, p, q, coeff, sum, EPS, epsilon, invsum
      PARAMETER (EPS = 1.E-8, epsilon = 1.E-30)
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL EXARET, EXARVO, EXANVA, PMV, PVX  
c     
c     intrinsic functions
      INTRINSIC SQRT       
c
c-----------------------------------------------------------------------
c
c     initialisation 
      info = 0
c      CALL IVX ( n, mvar, 0.0 )
      CALL IVX ( n, ivar, 0.0 )
      CALL IVX ( n, ipvar, 0.0 )
c
c     test confidence level
      IF ( (alpha .LT. EPS) .or. (alpha .GT. (1.-EPS)) ) THEN
         info = -102
         RETURN
      ENDIF
c
c     inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p  = alpha
      q  = 1. - p
      zc = dinvnr(p, q)
c
c     absolute portfolio variance
c     w' * Cov * w
      CALL OVTMCV ( n, cov, w, vol )
      IF (ABS(vol) .LT. epsilon) THEN
        DO i = 1,n
            ipvar(i) = w(i)
        ENDDO
        info = -1
        RETURN
      ENDIF
c
c     compute w'*cov
      CALL PMV ( n, n, cov, w, ivar )
c
c     compute  - Zc*(w'cov)/Volatility
      coeff = - zc/SQRT(vol)
      CALL PVX ( n, ivar, coeff, ivar )
c
c     marginal VaR
c     rho(i) - Zc*(w'cov)(i)/Volatility
      DO i = 1,n
        ivar(i) = rho(i) + ivar(i) 
      ENDDO
c
c     incremental contribution
      DO i = 1,n
        ivar(i) = w(i)*ivar(i)
      ENDDO
c
c     sum of incremental VaR (total VaR)
      sum = 0.
      DO i = 1, n 
        sum = sum + ivar(i)  
      ENDDO
      invsum = 1./sum    
c
c     incremental contribution in %  
      CALL PVX( n, ivar, invsum, ipvar )
c     
c     Normal VaR := max[0, min(1, -zc*volatility + mean)] 
      IF (sum .GT. 1.0) THEN
          CALL YV( n, ipvar, ivar )
      ENDIF
      IF (sum .LT. -1.0) THEN
        DO i =1,n
            ivar(i) = -ipvar(i)
        ENDDO
      ENDIF
c
      RETURN
      END
