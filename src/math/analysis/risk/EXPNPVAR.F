c=======================================================================
c
c     subroutine EXPNPVAR                                    
c
c     Ex-post Non Parametric Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPNPVAR ( n, m, x, alpha, dwork, npvar, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of return(s)                        integer
c            m      : number of non parametric point(s)          integer 
c            x      : return(s) (n)                               double
c            alpha  : probability level (0 < beta < 1)            double
c
c     WORKSPACE 
c            dwork  : vector (2*m)                                double
c
c     OUTPUT 
c            npvar  : non parametric Value-at-Risk                double
c            info   : = 0 successful exit                        integer
c
c     CALL 
c            EVMAX   : maximum of a vector
c            EVMIN   : minimum of a vector
c            KERNELD : Gaussian kernel method
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, m, info
      DOUBLE PRECISION alpha, npvar
      DOUBLE PRECISION x(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)      
c
c     local variables
      INTEGER i, p, pdy, pdz
      DOUBLE PRECISION eps, minx, maxx, sum       
      PARAMETER ( eps = 1.E-8 )
c
c     external subroutines
      EXTERNAL EVMAX, EVMIN, KERNELD
c     
c     intrinsic functions
      INTRINSIC min, max      
c
c-----------------------------------------------------------------------
c
c     initialization
      npvar = 0.0
      info  = 0
      p     = 0
      sum   = 0.0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdy = 1
c     pdy  : pointer , so ( m )
      pdz = pdy + m
c     pdz  : pointer , so ( m )
c
c     Total size of dwork array : 2*m
c
c----------------------------------------------------------------------
c
c     test confidence level: alpha
      IF ((alpha .LT. eps).or.(alpha .GT. (1.-eps))) THEN
         info = -102
         RETURN
      ENDIF
c
c     max/min of vector x
      CALL EVMAX ( n, x, maxx )
      CALL EVMIN ( n, x, minx )
c
c     Gaussian kernel method
      dwork(pdy) = minx
      DO i = 1, m-1
         dwork(pdy + i) = dwork(pdy + i - 1)  
     &                     + ((maxx - minx)/(m-1))
      ENDDO
      CALL KERNELD ( n, m, x, dwork(pdy), dwork(pdz), info )
      IF (info .LT. 0) RETURN
c     
c     non parametric Value-at-Risk at probability level alpha
      DO i = 1, m
          sum = sum + dwork(pdz + i - 1)
          IF (sum .LT. (1.0-alpha)) THEN
            p = p + 1   
          ENDIF
      ENDDO
      
      npvar = dwork(pdy + p)
c      npvar = min(npvar, 1.)
c      npvar = max(npvar, -1.)
c
      RETURN
      END

