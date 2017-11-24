c=======================================================================
c
c     subroutine CHECKFEASIT                                     
c
c     Checking non-emptiness of index tracking constraints
c
c-----------------------------------------------------------------------
      SUBROUTINE checkfeasit (n, w0, rho, rhob, neq, nin, ccst, bcst,
     &                        cinf, csup, delta, dwork, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                             integer
c            w0     : initial portfolio (n)                       double
c            rho    : mean returns (n)                            double
c            rhob   : mean return of index                        double
c            neq    : number of equality constraints             integer
c            nin    : number of inequality constraints           integer
c            ccst   : matrix of constraints (n*(neq+nin))         double
c            bcst   : vector of constraints (neq+nin)             double
c            cinf   : lower bound (n)                             double
c            csup   : upper bound (n)                             double
c            delta  : relative target return                      double
c
c     WORKSPACE 
c            dwork  : nin + neq                                   double
c                    
c     OUTPUT 
c            info   : diagnostic argument                        integer
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, info
      DOUBLE PRECISION w0(*), ccst(*), bcst(*), cinf(*), csup(*), 
     &                 rho(*), rhob, delta 
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables      
      INTEGER pdw
      DOUBLE PRECISION scal, EPS
      PARAMETER (EPS=1.E-8)
c
      pdw = 1
      info = 0     
      CALL checklinbox ( n, w0, neq, nin, ccst, bcst, cinf, csup,
     &                   dwork(pdw), info)
      CALL XV(n, w0, rho, scal)
      IF ((scal - rhob) + EPS .LT. delta) THEN
        info = -100
        RETURN    
      ENDIF 
      RETURN
      END
