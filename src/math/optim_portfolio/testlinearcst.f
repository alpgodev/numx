c=======================================================================
c
c     subroutine TESTLINEARCST 
c
c     This function tests if w(n) verifies the lower/upper bounds and
c     general linear constraints for a specified precision EPS=1.E-10
c
c     C*w <= b                      (general linear constraints) 
c     Cinf(i) <= w(i) <= Csup(i)    (lower/upper bounds)
c
c     if the constraints are verified then info = 0
c                                     else info = -115
c-----------------------------------------------------------------------
      SUBROUTINE testlinearcst ( n, w, neq, nin, ccst, bcst, cinf, csup,
     &                           dwork, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       w      : point (n)                                        double 
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       ccst   : matrix of constraints (nasset*(neq+nin))         double
c       bcst   : vector initial of constraints (neq+nin)          double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c
c     WORKSPACE 
c       dwork  : n                                                double
c
c     OUTPUT 
c       info   : diagnostic argument                             integer
c       
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION w(*), cinf(*), csup(*), ccst(n,*), bcst(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, m, pdy
      DOUBLE PRECISION a, b, EPS
      PARAMETER (EPS = 1.E-10)
c
c     external subroutines
      EXTERNAL PMTV
c
c-----------------------------------------------------------------------
c
c
c     initializations
      info = 0          ! diagnostic argument
      m    = neq + nin  ! number of general linear constraints
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdy = 1
c
c     Total size of dwork array = n
c
c-----------------------------------------------------------------------
c
c     test: cinf <= w <= csup
      DO i = 1,n
        IF ((w(i).LT.cinf(i)-EPS).OR.(w(i).GT.csup(i)+EPS)) THEN
            info = -115
            RETURN
        ENDIF
      ENDDO
c
c     no test todo if m=0
      IF (m .EQ. 0) THEN
        RETURN
      ENDIF
c
c     PMTV : matrix vector product, pdy(m) <- C'(n,m)*w(n) 
      CALL PMTV ( n, m, ccst, w, dwork(pdy) )  
c
c     test: C*w <= b
      DO i = 1,m      
        a = dwork(pdy+i-1)
        b = bcst(i)
        IF (a .GT. b+EPS) THEN
            info = -115
            RETURN
        ENDIF
      ENDDO
      RETURN
      END
