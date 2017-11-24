c=======================================================================
c
c     subroutine NAVEM                                        
c
c     Expectation-Maximization (EM) algorithm  
c     -> recover missing observations in time series 
c
c     NAV data structure 
c
c           n : number of points
c           p : number of assets
c
c     NAV(n,p) = | x(1,1), ..., x(1,p) | 
c                | x(2,1), ...,  miss  |
c                | x(3,1), ..., x(3,p) |
c                |   .   , ...,   .    |
c                |   .   , ...,   .    |
c                | x(n,1), ..., x(n,p) |
c
c     where x(2,p) = miss denote a missing NAV
c
c-----------------------------------------------------------------------
      SUBROUTINE navem ( n, p, NAV, miss, iwork, dwork, NAVc, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        n      : max. number of values (n > 1)                  integer
c        p      : number of assets (p > 0)                       integer
c        NAV    : NAV values (n*p)                                double
c        miss   : missing value (<= -1)                           double
c
c     WORKSPACE 
c        iwork  : 2*p                                            integer
c        dwork  : p*(13*p + 5*n + 10)                             double
c     
c     OUTPUT 
c        NAVc   : EM-recoved NAV values (n*p)                     double
c        info   : diagnostic argument                            integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION miss, NAV(n,*), NAVc(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, h, nh
      INTEGER piem, pdem, pdnav, pdR, pdRc
      DOUBLE PRECISION missret, r
c
c     external subroutines
      EXTERNAL YM, LOGRLACK, RETEM
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      CALL YV(n*p, NAV, NAVc)
      missret = -1000
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piem = 1
c     piem : pointer for RETEM (2*p)    
c
c     Total size of iwork array  = 2*p      
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdnav  = 1
c     pdnav  : pointer for NAV (n*p)
      pdR    = pdnav + ( n*p )
c     pdR    : pointer for log-returns (n*p)
      pdRc   = pdR + ( n*p )
c     pdRc   : pointer for corrected log-returns (n*p)
      pdem   = pdRc + ( n*p )
c     pdem   : pointer for RETEM (p*(13*p + 2*n + 10))
c
c     Total size of dwork array  = p*(13*p + 5*n + 10)
c             
c-----------------------------------------------------------------------
c
c     test if last date is not missing
      DO j = 1,p
        IF (NAV(n,j) .LE. miss) THEN
            info = -9004
            RETURN
        ENDIF
      ENDDO
c     
c     copy NAV -> NAVtmp
      CALL YV (n*p, NAV, dwork(pdnav))
      DO i = 1,n
        DO j = 1,p
            IF (NAV(i,j) .LE. miss) 
c             dwork(pdnav + p*(i-1) + j - 1) = missret
     &      dwork(pdnav + n*(j-1)+i-1) = missret
        ENDDO
      ENDDO
c
c     compute log-returns: log[NAV(t+h)/NAV(t)]
      h = 1
      nh = n - 1
      CALL logrlack ( n, p, dwork(pdnav), h, dwork(pdR), info )
c
c     EM-recover log-returns
      CALL retem ( nh, p, dwork(pdR), missret,
     &             iwork(piem), dwork(pdem), dwork(pdRc), info )
      IF (info .LT. 0) RETURN
c      
c     construct EM-recover NAV's
      DO i = 2,n
        DO j = 1,p
            IF (NAV(n - i + 1, j) .LE. miss) THEN
                r = dwork(pdRc + nh*(j-1) + nh - i + 1) ! returns
c                NAVc(n - i + 1,j) = NAVc(n - i + 2,j)/(1.0 + r)
                 NAVc(n - i + 1,j) = NAVc(n - i + 2,j)*EXP(-r)
            ENDIF
        ENDDO
      ENDDO
      RETURN
      END
