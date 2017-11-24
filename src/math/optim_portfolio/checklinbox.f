c=======================================================================
c
c     subroutine CHECKLINBOX                                   
c
c     Checking non-emptiness of linear constraints
c
c-----------------------------------------------------------------------
      SUBROUTINE checklinbox ( n, x0, neq, nin, ccst, bcst, cinf, csup,
     &                         dwork, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : problem size                                    integer
c       x0     : initial vector (n)                               double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : matrix of constraints (n*(neq+nin))              double
c       bcst   : vector of constraints (neq+nin)                  double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c
c     WORKSPACE 
c       dwork  : nin + neq                                        double
c                    
c     OUTPUT 
c       info   : diagnostic argument                             integer
c
c     CALL
c       PMTV
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, info
      DOUBLE PRECISION x0(*), cinf(*), csup(*), ccst(n,*), bcst(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     localvariables      
      INTEGER i, neqin, pdbt
      DOUBLE PRECISION temp, EPS 
      PARAMETER (EPS=1.E-8)
     
      pdbt = 1      
      info = 0      
c
      neqin = neq + nin
      CALL PMTV(n, neqin, ccst, x0, dwork(pdbt))
c      
      DO i = 1,neq
        temp = dwork(pdbt+i-1)-bcst(i)
        IF ((temp .LT. -EPS).OR.(temp .GT. EPS)) THEN
            info = -100
            RETURN
        ENDIF
      ENDDO
c      
      DO i = 1,nin
        IF (dwork(pdbt+neq+i-1).GT.(EPS+bcst(neq+i))) THEN
            info = -100
            RETURN
        ENDIF
      ENDDO
c      
      DO i = 1,n
        IF (x0(i)+EPS .LT. cinf(i)) THEN 
            info = -100
            RETURN
        ELSEIF (csup(i) .LT. x0(i)-EPS) THEN 
            info = -100
            RETURN
        ENDIF
      ENDDO  
      RETURN
      END
