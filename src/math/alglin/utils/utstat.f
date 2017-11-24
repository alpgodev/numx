c=======================================================================
c
c     Statistic Utilies                                      
c
c-----------------------------------------------------------------------
c
c        VARFAC  : explain variance by p factors
c        COVIND  : covariance of assets/index
c        COVR    : covariance matrix from centered returns
c        COVE    : exponentially weighted covariance matrix from centered returns
c        VARIAN  : variance
c        HVARIAN : variance on an history of assets
c        VOLAT   : volatility (standard deviation)
c        HVOLAT  : volatility on an history of assets
c        SEMVAR  : semi-variance ( variance(v), v(i) < mean(v(i)) )
c        SEMVOL  : semi-volatility
c        COVCOR  : conversion covariance  -> correlation
c        CORCOV  : conversion correlation -> covariance
c
c=======================================================================
c
c     subroutine VARFAC
c
c     explain variance by p factors
c
c-----------------------------------------------------------------------
      SUBROUTINE VARFAC ( N, eig, nfact, var, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : dimension of the eigenvalues vector             integer
c       eig(N) : sorted eigenvalues vector (N)                    double
c       nfact  : number of factors                               integer
c
c     OUTPUT 
c       var    : variance of nbfact factors                       double
c       info   : diagnostic argument                             integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c      
c     i/o parameters
      INTEGER N, nfact, info
      DOUBLE PRECISION var, eig(*)
c
c     local variables
      DOUBLE PRECISION vart, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     initialization
      info  = 0
c      
c     total variance (sum of the eigenvalues)
      CALL SEV ( N, eig, vart )
      IF ( ABS(vart) .LT. EPS ) THEN
         info = -1
         RETURN
      ENDIF
c
c     explained variance
      CALL SEV ( nfact, eig, var)
      var = var / vart
      RETURN
      END
c
c=======================================================================
c
c     subroutine COVIND
c
c     covariance between assets/index
c
c-----------------------------------------------------------------------
      SUBROUTINE COVIND ( p, n, X, Y, meanY, meanX, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p     : number of assets                                 integer
c       n     : number of returns                                integer
c       X     : assets matrix of returns (n*p)                    double
c       Y     : index vector of returns (n)                       double
c
c     OUTPUT 
c       meanY : mean return of index                              double
c       meanX : mean returns vector (n)                           double
c       cov   : covariance vector (n)                             double
c       info  : diagnostic argument                              integer
c
c     CALL   
c       MV    : mean of a vector
c       MCM   : mean of each column of a vectorized full matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER p, n, info
      DOUBLE PRECISION meanY, meanX(*), X(n,*), Y(*), cov(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION sum, tmp
c
c     initialization
      info = 0
c
c     test if ndate < 2      
      IF ( n .LT. 2 ) THEN
         info = -2
         RETURN
      ENDIF
c
c     mean of index and assets
      CALL MV ( n, Y, meanY )
      CALL MCM ( n, p, X, meanX )
c
c     cov = E(X)*E(Y) - E(X*Y)
      tmp = 1./FLOAT(n - 1)
      DO i = 1,p
         sum = 0.
         DO j = 1,n
            sum = sum + X(j,i) * Y(j)
         ENDDO
         cov(i) = tmp * ( sum - n*meanY*meanX(i) )
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine COVR
c
c     covariance matrix from centered returns 
c
c     cov(M,M) = (1/(N-1)) * X(N,M)'X(N,M) where X are centered returns
c
c-----------------------------------------------------------------------
      SUBROUTINE COVR ( N, M, X, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N        : number of rows                                integer
c       M        : number of columns                             integer
c       X(N*M)   : cebtered returns matrix (N*M)                  double
c
c     OUTPUT 
c       cov(M*M) : covariance matrix (M*M)                        double
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, M, info
      DOUBLE PRECISION alpha, X(N,*), cov(M,*)
c
c     initialization
      info = 0 
c      
c     test if N < 2  
      if ( N .LT. 2 ) THEN
         info = -2
         RETURN
      ENDIF
c
c     matrix product, cov(M,M) <- X'(N,M)*X(N,M)
      CALL PRMM ( N, M, X, cov )
c
c     scales a matrix by a constant, cov(M,M) <- alpha*cov(M,M)      
      alpha = 1./FLOAT(N-1)
      CALL PMX2 ( M, M, cov, alpha)
      RETURN
      END
c=======================================================================
c
c     subroutine COVE
c
c     exponential weighted covariance matrix from centered returns
c
c-----------------------------------------------------------------------
      SUBROUTINE COVE ( N, M, X, lambda, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N        : number of rows                                integer
c       M        : number of columns                             integer
c       X(N*M)   : centered returns matrix (N*M)                  double
c       lambda   : decay factor                                   double
c
c     OUTPUT 
c       cov(M*M) : exponential covariance matrix (M*M)            double
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, M, info
      DOUBLE PRECISION lambda, X(N,*), cov(M,*)
c
c     local variables
      INTEGER i, j, imat
      DOUBLE PRECISION sum, tmp, EPS
      PARAMETER( EPS = 1.E-4)
c
c     initialization
      info = 0 
c      
c     test N > 2
      IF ( N .LT. 2 ) THEN
         info = -2
         RETURN
      ENDIF
c
c     case lambda > 0.9999
      IF (lambda .GT. (1.-EPS)) THEN
         CALL COVR ( N, M, X, cov, info )
         RETURN
      ENDIF   
c     
      IF (lambda .LT. 1.E-15) THEN
         tmp = 1.0
      ELSE
        tmp = 0.0
        DO i = 1,N
           tmp = tmp + (lambda**(i-1))
        ENDDO
        tmp = 1./(tmp - 1.)
      ENDIF
      DO j = 1,M
         DO i = 1,j
            sum = 0.0
            DO imat = 1,N
               sum = sum + (lambda**(imat-1))*(X(imat,i)*X(imat,j))
            ENDDO
            cov(i,j) = tmp * sum
            IF ( i .NE. j ) cov(j,i) = cov(i,j)
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine VARIAN
c
c     Equally Weighted (Empirical) Variance
c
c-----------------------------------------------------------------------
      SUBROUTINE VARIAN ( N, X, var, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            N    : number of values                             integer
c            X(N) : input data, values (N)                        double
c
c     OUTPUT 
c            var  : variance                                      double
c            info : diagnostic argument                          integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, info
      DOUBLE PRECISION var, X(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION sum, sum2
c
c     initialization
      info = 0
      var  = 0.0
      sum  = 0.0
      sum2 = 0.0
c      
c     test if N < 2
      IF ( N .LT. 2 ) THEN
         info = -2
         RETURN
      ENDIF
c
c     square sum and sum of square 
      DO i = 1,N
         sum = sum + X(i)
         sum2 = sum2 + X(i)*X(i)
      ENDDO
c
c     variance
      var = ( sum2 - (sum*sum/N) ) / ( N - 1 )
      RETURN
      END
c
c=======================================================================
c
c     subroutine HVARIAN
c
c     Variance of P assets
c
c-----------------------------------------------------------------------
      SUBROUTINE HVARIAN ( N, P, X, var, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : number of values (N > 1)                        integer
c       P      : number of assets (P > 0)                        integer
c       X(N*P) : input data, values (N*P)                         double
c     OUTPUT 
c       var(P) : variance of each asset, vector (P)               double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       VARIAN : variance
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, P, info
      DOUBLE PRECISION var(*), X(N,*)
c
c     local variables
      INTEGER i
c
c     initialization
      info = 0
c      
      DO i = 1,P
         CALL VARIAN ( N, X(1,i), var(i), info )
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine VOLAT
c
c     Volatility (standard deviation)
c
c-----------------------------------------------------------------------
      SUBROUTINE VOLAT ( N, X, vol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : number of values (N > 1)                        integer
c       X(N)   : input data, values (N)                           double
c
c     OUTPUT 
c       vol    : volatility                                       double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       VARIAN : variance
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, info
      DOUBLE PRECISION vol, X(*)
c
c     local variables
      DOUBLE PRECISION var, EPS
      PARAMETER ( EPS = 1.E-30) 
c
c     initialization
      info = 0
      vol  = 0.0
c
c     variance
      CALL VARIAN ( N, X, var, info )
      IF (info .LT. 0) RETURN
c
c     test if variance is too small
      IF (var .LT. EPS) THEN
        vol = 0.0
        RETURN
      ELSE
c
c       volatility (standard deviation)
        vol = SQRT(var)
      ENDIF  
      RETURN
      END
c
c=======================================================================
c
c     subroutine HVOLAT
c
c     Volatility of P assets
c
c-----------------------------------------------------------------------
      SUBROUTINE HVOLAT ( N, P, X, vol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : number of dates (N > 1)                         integer
c       P      : number of assets (P > 0)                        integer
c       X(N*P) : input data, values (N*P)                         double
c
c     OUTPUT 
c       vol(P) : volatility of each asset, vector (P)             double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       VOLAT  : volatility
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, P, info
      DOUBLE PRECISION vol(*), X(N,*)
c
c     local variables
      INTEGER i
c
c     initialization
      info = 0
c
      DO i = 1,P
         CALL VOLAT ( N, X(1,i), vol(i), info )
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine SEMVAR
c
c     Semi-variance ( variance(v), for v(i) < mean(v(i)) )
c
c-----------------------------------------------------------------------
      SUBROUTINE SEMVAR ( N, X, svar, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N    : number of values (N > 1)                          integer
c       X(N) : input data, values (N)                             double
c
c     OUTPUT 
c       svar : semi-variance                                      double
c       info : diagnostic argument                               integer
c
c     CALL   
c       MV   : mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, info
      DOUBLE PRECISION svar, X(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION mu, sum, sum2, tmpx
c
c     initialization
      info = 0
      svar = 0.0
      sum  = 0.
      sum2 = 0.
      tmpx = 0.
c
c     test if N > 1      
      IF ( N .LT. 2. ) THEN
         info = -2
         RETURN
      ENDIF
c
c     mean of X 
      CALL MV ( N, X, mu )
c
      DO i = 1,N
         tmpx = MIN(X(i) - mu, 0.0)
         sum  = sum + tmpx
         sum2 = sum2 + (tmpx * tmpx)
      ENDDO
c
c     semi-variance of X
      svar = ( sum2 - (sum * sum / N) ) / ( N - 1. )
      RETURN
      END
c
c=======================================================================
c
c     subroutine SEMVOL
c
c     Semi-volatility ( volatility(v), v(i) < mean(v(i)) )
c
c-----------------------------------------------------------------------
      SUBROUTINE SEMVOL ( N, X, svol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N      : number of values (N > 1)                        integer
c       X(N)   : input data, values (N)                           double
c
c     OUTPUT 
c       svol   : semi-volatility                                  double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       SEMVAR : semi-variance ( variance(x), x(i) < mean(x(i)) )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, info
      DOUBLE PRECISION svol, X(*)
c
c     local variables
      DOUBLE PRECISION svar
c
c     initialization
      info = 0
      svol = 0.0
c
c     semi-variance
      CALL SEMVAR ( N, X, svar, info )
      IF (info .LT. 0) RETURN
c
c     semi-volatility (semi-standard deviation)
      svol = SQRT(svar)
      RETURN
      END
c
c=======================================================================
c
c     subroutine COVCOR
c
c     covariance matrix --> correlation matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE COVCOR ( N, cov, corr, vol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N         : size of the covariace matrix                 integer
c       cov(N*N)  : covariance matrix (N*N)                       double
c
c     OUTPUT 
c       corr(N*N) : correlation, matrix (N*N)                     double
c       vol(N)    : volatilities, vector (N)                      double
c       info      : diagnostic argument                          integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N, info
      DOUBLE PRECISION cov(N,*), corr(N,*), vol(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION vcov, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     initialization
      info = 0
c
c     extracting volatilities
      DO i = 1,N
         vcov = cov(i,i)
         IF ( vcov .GT. EPS ) THEN
            vol(i)= SQRT(vcov)
         ELSE
            info = -1
            RETURN
         ENDIF
      ENDDO
c
c     building the correlation matrix
      DO j = 1,N
         DO i = 1,j
            IF ( i .NE. j ) THEN
               corr(i,j) = cov(i,j)/( vol(i)*vol(j) )
               corr(j,i) = corr(i,j)
            ELSE
               corr(i,i) = 1.
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CORCOV
c
c     correlation matrix --> covariance matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE CORCOV ( N, corr, vol, cov )
c-----------------------------------------------------------------------
c
c     INPUT 
c       N         : size of the covariace matrix                 integer
c       corr(N*N) : correlation matrix (N*N)                      double
c       vol(N)    : volatilities vector (N)                       double
c
c     OUTPUT 
c       cov(N*N)  : covariance matrix (N*N)                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER N
      DOUBLE PRECISION corr(N,*), cov(N,*), vol(*)
c
c     local variables
      INTEGER i, j
c
      DO j = 1,N
         DO i = 1,j
            cov(i,j) = corr(i,j) * vol(i) * vol(j)
            IF ( i .NE. j ) cov(j,i) = cov(i,j)
         ENDDO
      ENDDO
      RETURN
      END
