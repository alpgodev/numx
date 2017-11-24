c=======================================================================
c
c     subroutine COVL                                        
c
c     Empirical Covariance matrix with different data length (cf. covm.f)
c
c-----------------------------------------------------------------------
      SUBROUTINE covl ( n, p, length, x, dwork, cov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values                           integer
c            p      : number of assets (p > 0)                   integer
c            length : lenght of each asset hist. (p)             integer
c            x      : assets values (n*p)                         double
c
c     WORKSPACE 
c            dwork  : (p*(1 + n + p) + 2*n)                       double
c     
c     OUTPUT 
c            cov    : covariance matrix (p*p)                     double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c           EVIMAX  : maximum of a vector (integers)
c           YMP     : copy a sub-block of a vectorized matrix in a vectorized matrix
c           COVMVM  : covariance matrix and the mean vector
c           VARIAN  : variance
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      INTEGER length(*)
      DOUBLE PRECISION x(n,*), cov(p,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, l, izero, iun, ideux, maxlength, lmin,
     &        pdwork, pdr, pdw, pdv, pdm
      PARAMETER ( izero = 0, iun = 1, ideux = 2 )
c     
c     intrinsic functions
      INTRINSIC min
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1
      pdr    = pdwork
c     pdr : pointer for sub-values, needs (n*2)
      pdw = pdr + ( 2*n )         
c     pdw : pointer for dwork COVMVM, needs ( n*p)
      pdv = pdw + ( n*p  )
c     pdv : pointer for rho (in COVM), needs (p)
      pdm = pdv + ( p )
c     pdm : pointer for covariance, needs (p*p)      
c
c     Total size of dwork array  = p*(1 + n + p) + 2*n
c             
c-----------------------------------------------------------------------
c
c     find the max length of assets history
      CALL EVIMAX ( p, length, maxlength )
      IF (maxlength .GT. n) THEN
         info = -100001
         RETURN
      ENDIF
c
c     covariance matrix
      DO i = 1,p-1
         DO j = i+1,p
c
c           covariance of two assets (on the smallest length)
            lmin  = min(length(i), length(j))
            CALL  YMP ( n, p, x, lmin, iun, iun, i,
     &                 dwork(pdr), info )
            CALL  YMP ( n, p, x, lmin, iun, iun, j,
     &                  dwork(pdr + lmin), info )
c           YMP : info = 0 by construction
            CALL covmvm ( lmin, ideux, dwork(pdr), dwork(pdw),
     &                  dwork(pdv), dwork(pdm), info)
            IF (info .lt. 0) RETURN
c
c           writing in full covariance matrix of non-diagonal elements
c           of 2x2 matrix
            cov(j,i) = dwork(pdm + 1)
            cov(i,j) = dwork(pdm + 2)
         ENDDO
      ENDDO
c
c     computing variance (with all values)
      DO i = 1,p
         l = length(i)
         CALL VARIAN ( l, x(1,i), cov(i,i), info )
         IF (info .LT. 0) RETURN
      ENDDO
      RETURN
      END
