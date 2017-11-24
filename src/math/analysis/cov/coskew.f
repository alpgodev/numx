c=======================================================================
c
c     subroutine COSKEW                     
c
c     Empirical CoSkewness matrix, p by (p*p)
c
c     Statistical measure that calculates the symmetry of a variable's 
c     probability distribution in relation to another variable's probability 
c     distribution symmetry. All else being equal, a positive coskewness 
c     means that the first variable's probability distribution is skewed 
c     to the right of the second variable's distribution.  
c  
c     In finance, coskewness can be used as a supplement to the covariance 
c     calculation of risk estimation. Usually, coskewness is calculated 
c     using a security's historic price data as the first variable, 
c     and the market's historic price data as the second. This provides 
c     an estimation of the security's risk in relation to market risk.
c
c     An investor would prefer a positive coskewness because this represents 
c     a higher probability of extreme positive returns in the security 
c     over market returns.
c         
c                           n
c                      1   ---            -                 -                 -
c     coskew(i,j,k) = ---   \   [X(t,i) - X(i)] * [X(t,j) - X(j)] * [X(t,k) - X(k)] 
c                      n    / 
c                           ---
c                           t=1
c     where
c     
c                  n
c       -      1  ---
c       X(j) = -  \   X(i,j) for j = 1,...,p
c              n  / 
c                 ---
c                 i=1
c
c     Note: In co-Skewness only (n+2)(n+1)n/6 elements must be computed.
c
c-----------------------------------------------------------------------
      SUBROUTINE coskew ( n, p, X, dwork, coskewness, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n           : number of value(s) (n > 1)                 integer
c       p           : number of asset(s) (p >= 1)                integer
c       X           : value(s) (n*p)                              double
c
c     WORKSPACE 
c       dwork       : vector( p*n +p)                             double
c
c     OUTPUT 
c        coskewness : empirical skewness matrix (p*p*p)           double
c        info       : diagnostic argument                        integer
c
c     CALL   
c        TENSORK    : estimation of co-Skewness[k] matrix
c        YM         : matrix copy
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(*), coskewness(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER k, pdtensor, pdHk
c
c     external subroutines
      EXTERNAL TENSORK, YM
c
c-----------------------------------------------------------------------
c
c     initialisations
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdtensor = 1
c     pdtensor: pointer for TENSORK (p)      
      pdHk     = pdtensor + ( p )
c     pdHk    : pointer for Hk tensor ( p*p )
c
c     Total size of dwork array = p*(p + 1) 
c
c----------------------------------------------------------------------
c
c     loop
      DO k = 1,p
c       tensor H[k]
        CALL TENSORK ( n, p, k, X, dwork(pdtensor), dwork(pdHk), info)
        IF (info .LT. 0) RETURN
c        
        CALL YM ( p, p, dwork(pdHk), coskewness(1+(k-1)*p*p) )        
      ENDDO
c
      RETURN
      END
