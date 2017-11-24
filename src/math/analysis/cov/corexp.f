c=======================================================================
c
c     subroutine COREXP                                        
c
c     Empirical Correlation matrix and mean vector
c
c                        cov(i,j)                              
c     corr(i,j) = -----------------------, 
c                 sqrt[cov(i,i)*cov(j,j)]
c
c     where cov is exponentially weighted covariance matrix
c
c                       n
c                       ---                              
c     cov(i,j) = (1-L)  \    L^(k-1)*[X(k,i) - rho(i)] * [X(k,j) - rho(j)] 
c                       / 
c                       ---
c                       k=1
c     
c                  n
c              1  ---
c     rho(j) = -  \   X(i,j) for j = 1,...,p
c              n  / 
c                 ---
c                 i=1
c
c     L : exponential parameter
c
c-----------------------------------------------------------------------
      SUBROUTINE corexp ( n, p, X, lambda, dwork, corr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values (n > 1)                   integer
c            p      : number of asset(s) (n >= 1)                integer
c            X      : value(s) (n*p)                              double
c            lambda : exponential parameter                       double
c
c     WORKSPACE 
c            dwork  : vector ( p*(n + p + 1) )                    double
c
c     OUTPUT 
c            corr   : exp. weighted correlation matrix (p*p)      double
c            info   : diagnostic argument                        integer
c 
c     CALL  
c            COVEXP : empiric exp. weighted covariance matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION lambda, X(n, *), corr(p,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdwork, pdw, pdcov, pdsqr, i, j, kii, kij
      DOUBLE PRECISION var, eps
      PARAMETER ( eps = 1.E-30 )
c
c     external subroutines
      EXTERNAL covexp
c     
c     intrinsic functions
      INTRINSIC SQRT
c
c-----------------------------------------------------------------------
c
c     initialization 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1
      pdw  = pdwork
c     pdw   : pointer for COVEXP, so ( p*(n+1) ) more
      pdcov = pdw + ( p*(n + 1) )
c     pdcov   : pointer for cov. matrix, so ( p*p ) more
      pdsqr = pdwork + ( p*p )
c     psqr   : pointer for square root of diagonal, so ( p ) more
c
c     Total size of dwork array = p*(n + p + 1) 
c
c-----------------------------------------------------------------------
c
c     empirical exponential weighted covariance matrix
      CALL covexp ( n, p, X, lambda, dwork(pdw), dwork(pdcov), info)
      IF (info .LT. 0) RETURN
c
c     standard deviation:  square root of diag(cov)
      DO i = 1,p
         kii = pdcov + (i - 1)*(p + 1)
         var = dwork(kii)
         IF ( var .GT. eps ) THEN
            dwork(pdsqr+i-1)= SQRT(var)
         ELSE
            info = -5
            RETURN
         ENDIF
      ENDDO
c
c     empirical correlation matrix
      corr(1,1) = 1.0
      DO j = 2,p
         DO i = 1,j
            kij       = pdcov + (i - 1)*p + (j - 1)
            corr(i,j) = dwork(kij)/( dwork(pdsqr+i-1)*dwork(pdsqr+j-1) )
            IF ( i .NE. j ) corr(j, i) = corr(i, j)
         ENDDO
         corr(j, j) = 1.0
      ENDDO
      RETURN
      END
