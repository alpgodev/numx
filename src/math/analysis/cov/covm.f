c=======================================================================
c
c     subroutine COVM                                      
c
c     Empirical Covariance matrix and mean vector
c
c                       n
c                 1    ---                              
c     cov(i,j) = ---   \   [X(k,i) - rho(i)] * [X(k,j) - rho(j)] 
c               (n-1)  / 
c                      ---
c                      k=1
c     
c                  n
c              1  ---
c     rho(j) = -  \   X(i,j) for j = 1,...,p
c              n  / 
c                 ---
c                 i=1
c
c-----------------------------------------------------------------------
      SUBROUTINE covm ( n, p, X, dwork, cov, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of value(s) (n > 1)                      integer
c       p      : number of asset(s) (p >= 1)                     integer
c       X      : value(s) (n*p)                                   double
c
c     WORKSPACE 
c       dwork  : p*n + p                                          double
c
c     OUTPUT 
c       cov    : empirical covariance matrix (p*p)                double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       MCM    : mean of each column of a vectorized full matrix
c       EMPCOV : empirical covariance matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(*), cov(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      integer pdwork, pdrho
c
c     external subroutines
      EXTERNAL MCM, empcov
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdrho = 1
c     prho: pointer ont the mean return rho, so p      
      pdwork = pdrho + p
c     pdwork : pointer for relatives rents (n*p)
c
c     Total size of dwork array = p*n + p 
c
c----------------------------------------------------------------------
c
c     means return(s)
      CALL MCM ( n, p, X, dwork(pdrho) )
c
c     empirical covariance matrix
      CALL empcov ( n, p, X, dwork(pdrho), dwork(pdwork), cov, info )
      RETURN
      END
