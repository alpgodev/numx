c=======================================================================
c
c     subroutine EMPCOV                                      
c
c     Empirical covariance matrix and mean:
c
c                         n
c                 1    ---                              
c     cov(i,j) = ---   \   [X(k,i) - rho(i)] * [X(k,j) - rho(j)] 
c               (n-1)  / 
c                      ---
c                      k=1
c
c-----------------------------------------------------------------------
      SUBROUTINE empcov ( n, p, X, rho, XR, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of value(s) (n > 1)                       integer
c       p     : number of asset(s) (p >= 1)                      integer
c       X     : matrix of values (n*p)                            double
c       rho   : mean vector (p)                                   double
c
c     OUTPUT 
c       XR    : relative [X(.,j) - rho(j)] (n*p)                  double
c       cov   : empirical covariance matrix (p*p)                 double
c       info  : diagnostic argument                              integer
c
c     CALL   
c        COVR : covariance matrix from centred values : X - E(X)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(n,*), rho(*), cov(p,*), XR(n,*)
c
c     local variables 
      INTEGER i, j
      DOUBLE PRECISION mean
c
c     external subroutines
      EXTERNAL covr
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0     
c      
c     relative return(s)
      DO j = 1,p
         mean = rho(j)
         DO i = 1,n
            XR(i, j) = X(i,j) - mean
         ENDDO
      ENDDO
c
c     cov(p,p) <- [XR'(n,p)*XR(n,p)]/(n-1)
      CALL covr ( n, p, XR, cov, info )
      RETURN
      END
