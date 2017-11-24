c=======================================================================
c
c     subroutine COKURT
c
c     Empirical CoKurtosis matrix, p by (p*p)
c         
c                            n
c                        1  ---         -             -             -             -
c     cokurt(i,j,k,l) = --- \   [X(t,i)-X(i)]*[X(t,j)-X(j)]*[X(t,k)-X(k)]*[X(t,l)-X(l)]
c                        n  / 
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
c-----------------------------------------------------------------------
      SUBROUTINE cokurt ( n, p, X, dwork, cokurtosis, info)
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
c        cokurtosis : empirical cokurtosis matrix (p*p*p*p)       double
c        info       : diagnostic argument                        integer
c
c     CALL   
c        TENSORK    : estimation of co-Skewness[k] matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(*), cokurtosis(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER k, l, pdtensor, pdHk
c
c     external subroutines
      EXTERNAL TENSORKK, YM
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdtensor = 1
c     pdtensor: pointer for TENSORKK (p)      
      pdHk     = pdtensor + ( p )
c     pdHk    : pointer for Hk tensor ( p*p )
c
c     Total size of dwork array = p*(p + 1) 
c
c----------------------------------------------------------------------
c
c     loop
      DO k = 1,p
        DO l = 1,p
c           tensor H[k]
            CALL TENSORKK (n,p,k,l,X,dwork(pdtensor),dwork(pdHk),info)
            IF (info .LT. 0) RETURN
c        
            CALL YM ( p, p, dwork(pdHk), 
     &                cokurtosis(1+(l-1)*p*p + (k-1)*p*p*p) )      
        ENDDO
      ENDDO
c
      RETURN
      END
