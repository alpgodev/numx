c=======================================================================
c
c     subroutine COVLACK                                        
c
c     Empirical Covariance matrix with missing data (lack of data)
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
      SUBROUTINE covlack ( n, p, x, z, dwork, cov, ind, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : max. number of values                      integer
c            p      : number of assets (p > 1)                   integer
c            x      : asset(s) return(s) values (n*p)             double
c            z      : backup cov. matrix (p*p)                    double
c
c     WORKSPACE 
c            dwork  : n(p + 4) + 6                                double
c     
c     OUTPUT 
c            cov    : covariance matrix (p*p)                     double
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
      DOUBLE PRECISION x(n,*), z(p,*), cov(p,*)
      INTEGER ind(p,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, k, m, l, nb, pdwork, pdcovm, pdcov, pdrho, pdx1, 
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
      pdcovm = pdwork
c     pdcovm : pointer COVMVM workspace (n*p)
      pdcov  = pdcovm + ( n*p )
c     pdcov  : pointer for tmp cov (4)       
      pdrho  = pdcov + ( 4 ) 
c     pdrho  : pointer for rho vector (COVM), (2)
      pdx1   = pdrho + ( 2 )
c     pdx1   : pointer for X1, max(n)
      pdx2   = pdx1 + ( n )
c     pdx2   : pointer for X2, max(n)
      pdx    = pdx2 + ( n ) 
c     pdx    : pointer for X = [X1, X2], max( 2*n )       
c
c     Total size of dwork array  = n(p + 4) + 6 +p
c             
c-----------------------------------------------------------------------
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
c           case m=0,1, no enough data for computation           
            IF (m .LT. 2) THEN
                cov(k,j) = z(k,j)
                cov(j,k) = cov(k,j)
                ind(k,j) = 1
                ind(j,k) = 1
                out = out + 1    
            ELSE
                ind(k,j) = 0
                ind(j,k) = 0
c
c               X = [x1, x2] with X(m,2)
                DO i = 1,m
                    dwork(pdx + i - 1)      = dwork(pdx1 + i - 1) 
                    dwork(pdx + m + i - 1)  = dwork(pdx2 + i - 1) 
                ENDDO
c
                CALL covmvm ( m, nb, dwork(pdx), dwork(pdcovm),
     &                     dwork(pdrho), dwork(pdcov), info)
                IF (info .LT. 0) RETURN
                cov(k,j) = dwork(pdcov + 1)
                cov(j,k) = cov(k,j)
             ENDIF
c        
        ENDDO
      ENDDO        
c     
c     computing variance cov(i,i)
      DO j = 1,p
        l = 0
        DO i = 1,n
            IF (x(i,j) .GT. missing) THEN
                l = l + 1
                dwork(pdx1 + l - 1) = x(i,j)
            ENDIF
        ENDDO
c
c       case l=0,1  no enough data for computation
        IF (l .LT. 2) THEN
            cov(j,j) = z(j,j)
            ind(j,j) = 1
            out = out + 1
        ELSE
            ind(j,j) = 0
            CALL VARIAN ( l, dwork(pdx1), cov(j,j), info )
        ENDIF
        IF (info .LT. 0) RETURN
      ENDDO
      RETURN
      END
