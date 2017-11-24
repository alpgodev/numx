c=======================================================================
c
c     subroutine VOLCONTRIB                                  
c
c     Volatility Contribution (marginal risk contribution)
c
c     vol(i) = w(i)*[w'*Cov] / (Total Volatility)
c            = w(i)*[w'*Cov] / sqrt(w'*cov*w)
c
c-----------------------------------------------------------------------
      SUBROUTINE VOLCONTRIB( n, cov, w, mvol, mper, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c       w      : weight(s) (n)                                    double 
c       cov    : covariance matrix (n*n)                          double
c
c     OUTPUT 
c       mvol   : marginal volatility contribution (n)             double
c       mper   : marginal volatility contribution in % (n)        double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       IVX
c       OVTMCV    
c       PMV    : Matrix * Vector = Vector
c       PVX    : Vector * Scalar = Vector
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION w(*), cov(*), mvol(*), mper(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION var, sum, EPS
      PARAMETER (EPS = 1.E-30)
c
c     external subroutines
      EXTERNAL OVTMCV, PMV, PVX  
c
c     intrinsic
      INTRINSIC SQRT      
c
c-----------------------------------------------------------------------
c
c     initialisation 
      info = 0
      CALL IVX ( n, mvol, 0.0 )
      CALL IVX ( n, mper, 0.0 )     
c
c     absolute portfolio variance
c     w' * Cov * w
      CALL OVTMCV ( n, cov, w, var )
c
c     w'*cov
      CALL PMV ( n, n, cov, w, mvol )
c
c     marginal contribution w(i)*[(w'cov)/volatility]
      DO i = 1,n
        mvol(i) = w(i)*mvol(i)/SQRT(var)
      ENDDO
      sum = 0.
      DO i = 1,n
        sum = sum + mvol(i)
      ENDDO
      sum = 1./sum
c
c     marginal contribution in % (marginal contribution / volatility)  
      CALL PVX ( n, mvol, sum, mper )
c
      RETURN
      END  
