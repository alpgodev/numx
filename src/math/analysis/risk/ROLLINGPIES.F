c=======================================================================
c
c     subroutine ROLLINGPIES                               
c
c     Historical weights (rolling pies).
c
c     w[i] = w[0] * (x[i]/x[0])*(y[0]/y[i])
c
c     w[t] : weight at time t
c     x[t] : asset value at time t
c     y[t] : portfolio value at time t
c
c-----------------------------------------------------------------------
      SUBROUTINE ROLLINGPIES ( n, p, wini, x, y, mode, dwork, w, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of value(s) (n > 1)                       integer
c       p     : number of asset(s) (p > 0)                       integer 
c       wini  : initial weight (p)                                double
c       x     : asset value(s) (n*p)                              double
c       y     : portfolio value(s) (n)                            double
c       mode  : =1 past values, =2 future values                 integer
c
c     WORKSPACE 
c       dwork : 2*n                                               double
c
c     OUTPUT 
c       w     : historical weight(s) (n*p)                        double
c       info  : = 0 successful exit                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, mode, info
      DOUBLE PRECISION wini(*), x(n,*), y(*), w(n,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables      
      INTEGER i, j, pdx, pdw
      DOUBLE PRECISION wi, sum, EPS
      PARAMETER (EPS = 1.E-30)
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0  
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdx = 1
c     pdx : pointer for temporary X, needs (n)              
      pdw = pdx + ( n )
c     pdwo  : pointer for temporary W, need (n)
c
c     Total size of dwork array = 2*n
c
c-----------------------------------------------------------------------  
c
c     test number of data/assets
      IF (n .LE. 1) THEN
        info = -2
        RETURN 
      ENDIF
      IF (p .LT. 1) THEN
        RETURN 
      ENDIF
c
c     loop on each asset
      DO i = 1,p
c
c       initial weight of asset i      
        wi = wini(i)
c
c       copy the i-column of X(n,p) in V(n)         
        CALL YCMV ( n, p, x, i, dwork(pdx), info )
c
c       compute the vector of weights for asset i
        CALL ROLPIES ( n, wi, dwork(pdx), y, mode, dwork(pdw), info)
        IF (info .LT. 0) RETURN
c
c       copy the vector of weights in matrix W
        CALL YVCM ( n, p, w, dwork(pdw), i, info )            
      ENDDO
c
c     weights normalisation
      IF (p .GT. 1) THEN
      DO i = 1,n
        sum = 0.0
c        
        DO j = 1,p
            sum = sum + w(i,j)
        ENDDO
c        
        IF (sum .LT. EPS) THEN
         info = -9001
         RETURN
        ENDIF
c        
        DO j = 1,p
            w(i,j) = w(i,j)/sum
        ENDDO
      ENDDO  
      ENDIF
c
      RETURN
      END
