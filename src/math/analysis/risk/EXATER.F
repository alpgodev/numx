c=======================================================================
c
c     subroutine EXATER                                      
c
c     Ex-ante Tracking Error
c
c-----------------------------------------------------------------------
      SUBROUTINE EXATER (n, cov, covben, varben, w, dwork, tracke, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                             integer
c       w      : weight(s) (n)                                   double
c       cov    : covariance matrix (n*n)                         double
c       covben : covariance of benchmark (n)                     double
c       varben : variance of benchmark                           double
c
c     WORKSPACE 
c       dwork  : n*(2*n + 3) + 2                                 double
c
c     OUTPUT 
c       tracke : ex-ante tracking error                          double
c       info   : diagnostic argument                            integer
c                    
c     CALL   
c       OVTMCV : V'*M*V = scalar
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION varben, tracke, w(*), cov(*), covben(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      integer pcov, px, pmatq, i, j, kc, kq 
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL OVTMCV
c     
c     intrinsic functions
      INTRINSIC SQRT
c
c-----------------------------------------------------------------------
c
c     initialization
      tracke = 0.0
      info   = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pcov = 1
c     pcov  : pointer for covariance matrix (n*n)
      px   = pcov + ( n*n )
c     px    : pointer for vextor x, so (n+1) more
      pmatq = px + ( n+1 )
c     pmatq : pointer for matrix Q, so (n+1)*(n+1) more
c
c     Total size of dwork array = n*(2*n + 3) + 2
c
c-----------------------------------------------------------------------
c
      CALL YM ( n, n, cov, dwork(pcov) )
c
c     construction of vector X = [ weight , -1 ]
      DO i = 1,n
         dwork(px+i-1) = w(i)
      ENDDO
      dwork(px+n) = -1
c
c     construction of matrix Q = [ cov     covben' ]
c                                [ covben  varben  ]
      kc = pcov
      kq = pmatq
      DO i = 1,n
         DO j = 1,n
            dwork(kq) = dwork(kc)
            kc = kc + 1
            kq = kq + 1
         ENDDO
         dwork(kq) = covben(i)
         kq = kq + 1
      ENDDO
      DO i = 1,n
         dwork(kq) = covben(i)
         kq = kq + 1
      ENDDO
      dwork(kq) = varben
c
c     tracking error = sqrt(X'*Q*X)
      CALL OVTMCV ( n+1, dwork(pmatq), dwork(px), tracke )
      IF (tracke .LT. EPS) THEN
        tracke = 0.0
        RETURN
      ENDIF  
      tracke = SQRT(tracke)
c
      RETURN
      END
