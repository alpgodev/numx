c=======================================================================
c
c     subroutine EXPNPCVAR                                   
c
c     Ex-post Non Parametric Conditional Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPNPCVAR ( n, m, x, alpha, dwork, npcvar, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of return(s)                             integer
c       m      : number of non parametric point(s)               integer 
c       x      : return(s) (n)                                    double
c       alpha  : probability level (0 < alpha < 1)                double
c
c     WORKSPACE 
c       dwork  : vector (2*m)                                     double
c
c     OUTPUT 
c       npcvar : non parametric Conditional Value-at-Risk         double
c       info   : diagnostic argument                             integer
c
c     CALL 
c      EVMAX   : maximum of a vector
c      EVMIN   : minimum of a vector
c      KERNELD : Gaussian kernel method
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, m, info
      DOUBLE PRECISION alpha, npcvar
      DOUBLE PRECISION x(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)      
c
c     local variables
      INTEGER i, p, pdy, pdz
      DOUBLE PRECISION EPS, minx, maxx, npvar, sum, tmp       
      PARAMETER ( EPS = 1.E-8 )
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
      info   = 0
      npvar  = 0.0
      npcvar = 0.0
      p      = 0
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
      IF ((alpha .LT. EPS).or.(alpha .GT. (1.-EPS))) THEN
         info = -102
         RETURN
      ENDIF
c
c     max/min of vector x
      CALL EVMAX ( n, x, maxx )
      CALL EVMIN ( n, x, minx )
      minx = minx - (maxx - minx)*0.20
c
c     Gaussian non parametric kernel method 
      dwork(pdy) = minx
      DO i = 1, m-1
         dwork(pdy + i) = dwork(pdy + i - 1)
     &                     + ((maxx - minx)/(m-1))
      ENDDO
      CALL KERNELD ( n, m, x, dwork(pdy), dwork(pdz), info )
      IF (info .LT. 0) RETURN
c
c     non parametric Conditionnal Value-at-Risk at probability level alpha%
      sum = 0.0
      DO i = 1, m
          sum = sum + dwork(pdz + i - 1)
      ENDDO
      tmp = 0.0
      DO i = 1, m
          tmp = tmp + dwork(pdz + i - 1)/sum
          IF (tmp .LT. ABS(1.0 - alpha)) THEN
            p = p + 1
          ENDIF  
      ENDDO
      npvar = dwork(pdy + p)
c
c     non parametric CVaR at probability level alpha%
      sum = 0.0
      IF (p .GT. 0) THEN
        DO i = 1,p
            sum = sum + dwork(pdy + i - 1)
        ENDDO
        npcvar = sum / p
      ELSE
        info = -2  
      ENDIF
c
      RETURN
      END
