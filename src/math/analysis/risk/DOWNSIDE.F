c=======================================================================
c
c     subroutine DOWNSIDE                                    
c
c     Downside Risk
c
c     downside := E[min(x(i) - mu, 0)]^k for k = 1,2,3,4
c
c       k = 1: semi-mean E[min(x(i) - mu, 0)]
c       k = 2: semi-variance E[min(x(i) - mu, 0)]^2         
c              semi-volatility sqrt{E[min(x(i) - mu, 0)]^2} 
c       k = 3: semi-skewness
c       k = 4: semi-kurtosis
c
c-----------------------------------------------------------------------
      SUBROUTINE DOWNSIDE ( n, k, mu, x, dwork, dwnside, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : number of values ( n > 1)                  integer
c            k      : k-moment (k=1,2,3,4)                       integer
c            mu     : threshold return                            double
c            x      : values (n)                                  double
c
c     WORKSPACE :
c            dwork : (n)                                         integer 
c    
c     OUTPUT :
c            dwnside  : downside risk                             double
c            info     : diagnostic argument                      integer
c                       = 0 computation performed
c                       if -2, k is not equal to 1,2,3,4  
c    
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n, k, info
      DOUBLE PRECISION mu, dwnside, x(*)
c
c     workspace
      DOUBLE PRECISION dwork(*)      
c
c     local variables
      INTEGER i, pdwork
      DOUBLE PRECISION meany, vary 
c
c-----------------------------------------------------------------------
c
c     initialization
      info     = 0
      dwnside  = 0.0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwork = 1
c     pdwork  : pointer, so (n)      
c
c----------------------------------------------------------------------
c
c     test k
      IF ((k.le.0).or.(k.ge.5)) THEN
        info = -2
        RETURN
      ENDIF
c
c     y := min[x(i) - mu, 0]
      DO i = 1,n 
         dwork(pdwork + i - 1) = min(x(i) - mu, 0.0)
      ENDDO
c
c     case k = 1, semi-mean (with threshold mu)
      IF (k .eq. 1) THEN
        CALL MV ( n, dwork(pdwork), dwnside )
      ENDIF  
c
c     case k = 2, semi-variance (with threshold mu)
      IF (k .eq. 2) THEN      
        CALL VARIAN ( n, dwork(pdwork), dwnside, info )
      ENDIF      
c
c     case  k = 3, semi-skewness (with threshold mu)
      IF (k .eq. 3) THEN      
        CALL MV ( n, dwork(pdwork), meany )
        CALL VARIAN ( n, dwork(pdwork), vary, info ) 
        CALL SKEWN(n, dwork(pdwork), meany, vary, dwnside, info)
      ENDIF
c     
c     case k = 4, semi-kurtosis (with threshold mu) 
      IF ( k .eq. 4) THEN
        CALL MV ( n, dwork(pdwork), meany )
        CALL VARIAN ( n, dwork(pdwork), vary, info ) 
        CALL KURTOS(n, dwork(pdwork), meany, vary, dwnside, info)
      ENDIF    
c
      RETURN
      END
