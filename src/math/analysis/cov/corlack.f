c=======================================================================
c
c     subroutine CORLACK                                        
c
c     Empirical correlation matrix with missing data (lack of data)
c
c     Returns data structure
c
c     X(n,p) = | x(1,1), ..., x(1,p) |  
c              | x(2,1), ..., -1000  |
c              | x(3,1), ..., x(3,p) |
c              |   .   , ...,   .    |
c              |   .   , ...,   .    |
c              | x(n,1), ..., x(n,p) |
c
c     where x(2,p) = -1000 denote a missing return
c
c-----------------------------------------------------------------------
      SUBROUTINE corlack ( n, p, x, z, dwork, corr, ind, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : max. number of values                      integer
c            p      : number of assets (p > 1)                   integer
c            x      : asset(s) return(s) values (n*p)             double
c            z      : backup corr. matrix (p*p)                   double
c
c     WORKSPACE 
c            dwork  : p*(n+p+1) + 4*n + 6                         double
c     
c     OUTPUT 
c            corr   : correlation matrix (p*p)                    double
c            ind    : index of used backup value(s) (p*p)        integer
c            info   : diagnostic argument                        integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION x(n,*), z(p,*), corr(p,*)
      INTEGER ind(p,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, k, m, nb, pdwork, pdcorm, pdcorr, pdrho, pdx1, 
     &        pdx2, pdx, out
      DOUBLE PRECISION missing
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      out = 0
      missing = -1000.0
      nb = 2
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1
      pdcorm = pdwork
c     pdcorm : pointer CORM workspace (p*(n + p + 1)
      pdcorr = pdcorm + ( p*(n + p + 1) )
c     pdcorr : pointer for tmp corr (4)       
      pdrho  = pdcorr + ( 4 ) 
c     pdrho  : pointer for rho vector (COVM), (2)
      pdx1   = pdrho + ( 2 )
c     pdx1   : pointer for X1, max(n)
      pdx2   = pdx1 + ( n )
c     pdx2   : pointer for X2, max(n)
      pdx    = pdx2 + ( n ) 
c     pdx    : pointer for X = [X1, X2], max( 2*n )       
c
c     Total size of dwork array  = p(n+p+1) + 4*n + 6
c             
c-----------------------------------------------------------------------
c     
c     correlation diagonal
      DO j = 1,p
        corr(j,j) = 1.     
      ENDDO 
c
c     loop on j=1,...,p-1 and k=j+1,...,p 
      DO j = 1,(p-1)
        DO k = (j+1),p
c
c           covariance of two assets (without missing returns)
            m = 0
            DO i = 1,n
                IF ((x(i,j) .GT. missing) .AND. (x(i,k) .GT. missing)) 
     &          THEN
                    m = m + 1
                    dwork(pdx1 + m - 1) = x(i,j) 
                    dwork(pdx2 + m - 1) = x(i,k)
                ENDIF
            ENDDO
c
c           case m=0,1 no enough data for computation           
            IF (m .LT. 2) THEN
                corr(k,j) = z(k,j)
                corr(j,k) = corr(k,j) 
                ind(k,j) = 1
                ind(j,k) = 1
                out = out + 1  
            ELSE
c
c               X = [x1, x2] with X(m,2)
                DO i = 1,m
                    dwork(pdx + i - 1)      = dwork(pdx1 + i - 1)
                    dwork(pdx + m + i - 1)  = dwork(pdx2 + i - 1)
                ENDDO
c             
c               covariance (i,j)
                CALL cormvm ( m, nb, dwork(pdx), dwork(pdcorm),
     &                      dwork(pdrho), dwork(pdcorr), info)
                IF (info .LT. 0) RETURN
                corr(k,j) = dwork(pdcorr + 1)
                corr(j,k) = corr(k,j)
            ENDIF
        ENDDO
      ENDDO       
c      info = out       
c      
      RETURN
      END
