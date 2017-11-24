c=======================================================================
c
c     subroutine COVEXP                                      
c
c     (Empirical) Exponentially Weighted Covariance matrix
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
      SUBROUTINE covexp ( n, p, X, lambda, dwork, cov, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of value(s) (n > 1)                 integer
c            p      : number of asset(s) (p >= 1)                integer
c            X      : value(s) (n*p)                              double
c            lambda : exponential parameter                       double
c
c     WORKSPACE 
c            dwork  : vector( p*(n + 1) )                         double
c
c     OUTPUT 
c            cov    : exp. weighted covariance matrix (p*p)       double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            MCM    : mean of each column of a vectorized full matrix
c            EMPCOV : empirical covariance matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION lambda, X(n, *), cov(p, *)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      integer pdwork, pdrho
c
c     external subroutines
      EXTERNAL MCM, empcovexp
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1
c     pdwork : pointer for relatives rents, so ( n*p ) more
      pdrho = pdwork + ( n*p )
c     pdrho  : pointer for mean return rho (p)
c
c     Total size of dwork array = p*(n + 1) 
c
c----------------------------------------------------------------------
c
c     means return(s)
      CALL MCM ( n, p, X, dwork(pdrho) )
c
c     empirical exp. weighted covariance matrix
      CALL empcovexp ( n, p, lambda, X, dwork(pdrho), dwork(pdwork),
     &                 cov, info )
c
      RETURN
      END

