c=======================================================================
c
c     subroutine EMPCOVEXP                                   
c
c     Empirical Exponentially Weighted Covariance matrix and mean
c
c                         n
c                        ---                              
c     cov(i,j) = (1 - L) \   L^(k-1)*[X(k,i) - rho(i)] * [X(k,j) - rho(j)] 
c                        / 
c                        ---
c                        k=1
c
c-----------------------------------------------------------------------
      SUBROUTINE empcovexp ( n, p, lambda, X, rho, XR, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n     : number of value(s) (n > 1)                  integer
c            p     : number of asset(s) (p >= 1)                 integer
c            lambda: exponetial parameter                        integer
c            X     : matrix of values (n*p)                       double
c            rho   : mean vector (n)                              double
c
c     OUTPUT 
c            XR    : relative [X(.,j) - rho(j)] (n*p)            double
c            cov   : empirical covariance matrix (p*p)           double
c            info  : diagnostic argument                        integer
c
c     CALL   
c        cove      : covariance matrix from centred values : X - E(X)
c                   (input matrix n*m, result : full symetric matrix m*m)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION lambda, X(n,*), rho(*), cov(p,*), XR(n,*)
c
c     local variables 
      INTEGER i, j
      DOUBLE PRECISION mean
c
c     external subroutines
      EXTERNAL cove
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
c     computation of sum on k of XR(k,i) * XR(k,j)
      CALL cove ( n, p, XR, lambda, cov, info )
c
      RETURN
      END

