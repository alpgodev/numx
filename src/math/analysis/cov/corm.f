c=======================================================================
c
c     subroutine CORM                                        
c
c     Empirical Correlation matrix and mean vector:
c
c                        cov(i,j)                              
c     corr(i,j) = -----------------------, cov is the empirical covariance
c                 sqrt[cov(i,i)*cov(j,j)]
c
c                  n
c              1  ---
c     rho(j) = -  \   X(i,j) for j = 1,...,p
c              n  / 
c                 ---
c                 i=1
c
c-----------------------------------------------------------------------
      SUBROUTINE corm ( n, p, X, dwork, corr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values (n > 1)                   integer
c            p      : number of asset(s) (n >= 1)                integer
c            X      : value(s) (n*p)                              double
c
c     WORKSPACE 
c            dwork  : vector ( p*(n + p + 2) )                    double
c
c     OUTPUT 
c            corr   : empirical correlation matrix (p*p)          double
c            info   : diagnostic argument                        integer
c 
c     CALL  
c            MCM    : mean of each column of a vectorized full matrix
c            EMPCOV : empiric covariance matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(n, *), corr(p,*)
c
c     workspaces     
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdwork, pcov, psqr, i, j, kij, kii, pdrho
      DOUBLE PRECISION val, myzero
      PARAMETER ( myzero = 1.E-30 )
c
c     external subroutines
      EXTERNAL MCM, empcov
c     
c     intrinsic functions
      INTRINSIC sqrt            
c
c-----------------------------------------------------------------------
c
c     initialization 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdrho = 1
c     prho: pointer ont the mean return rho, so p      
      pdwork = pdrho + p
c     pdwork : pointer for relatives rents, so ( n*p ) more
      pcov = pdwork + ( n*p)
c     pcov   : pointer for vectorized matrix, so ( p*p ) more
      psqr = pcov + ( p*p )
c     psqr   : pointer for square root of diagonal, so ( p ) more
c     pdnext = psqr + ( p )
c     pdnext: pointer of the next dwork array
c
c     Total size of dwork array = p*(n + p + 1) 
c
c-----------------------------------------------------------------------
c
c     mean returns
      CALL MCM ( n, p, X, dwork(pdrho))
c
c     empirical covariance matrix
      CALL empcov ( n, p, X, dwork(pdrho), dwork(pdwork), dwork(pcov),
     &  info )
      IF (info .LT. 0) RETURN
c
c     standard deviation:  square root of Diag(cov)
      DO i = 1,p
         kii = pcov + (i - 1)*(p + 1)
         val = dwork(kii)
         IF ( val .gt. myzero ) THEN
            dwork(psqr+i-1)= sqrt(val)
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
            kij       = pcov + (i - 1)*p + (j - 1)
            corr(i,j) = dwork(kij)/( dwork(psqr+i-1)*dwork(psqr+j-1) )
            IF ( i .ne. j ) corr(j, i) = corr(i, j)
         ENDDO
         corr(j,j) = 1.0
      ENDDO
c
      RETURN
      END

