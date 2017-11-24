c=======================================================================
c
c     subroutine ALLOCIT                                     
c
c     Index Tracking allocation srtategy
c
c     Min[ Tracking Error ]
c      s.t. 
c     rho*w - index perf. >= delta      (relative target return)
c     C*w <= b                          (linear constraints) 
c     Cinf <= w <= Csup                 (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c----------------------------------------------------------------------
      SUBROUTINE allocit ( n, cov, rho, covb, rhob,
     &                     neq, nin, ccst, bcst, cinf, csup, delta,
     &                     iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset                                 integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : expected returns (n)                             double
c       covb   : covariance assets-index (n)                      double
c       rhob   : index expected return                            double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : matrix of constraints (n*(neq+nin))              double
c       bcst   : vector of constraints (neq+nin)                  double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       delta  : relative target return                           double
c
c     WORKSPACE
c       iwork  : 3*n+2*nin+neq+6                                 integer 
c       dwork  : n*(2*n+2*nin+2*neq+18)+4*neq+6*nin+15            double
c                    
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       OPIT    : computing optimization for ALLOCIT (cf. OPIT.F)
c       TESTSDP : test if covariance matrix is SDP
c
c-----------------------------------------------------------------------   
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION rhob, delta
      DOUBLE PRECISION cov(*), rho(*), covb(*), cinf(*), csup(*), 
     &                 wopt(*), ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pisdp, piwo, pdsdp, pdwo
c
c     external subroutines
      EXTERNAL testsdp, opit
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for workspace of TESTSDP (n)
      piwo  = 1
c     piwo  : pointer for OPIT (3*n+2*nin+neq+6)
c
c     Total size of iwork array = 3*n + 2*nin + neq + 6 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))      
      pdwo  = 1
c     pdwo : pointer for OPIT
c            n*(2*n + 2*nin + 2*neq + 18) + 4*neq + 6*nin + 15
c
c     Total size of dwork array = n*(2*n+2*nin+2*neq+18)+4*neq+6*nin+15
c
c--------------------------------------------------------------------
c
c     covariance matrix SDP ?
      CALL testsdp (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF   
c
c     optimization - traking error minimazion
      CALL opit ( n, cov, covb, rho, rhob, delta,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info)
      IF (info .ge. 0) THEN
        CALL checkfeasit (n, wopt, rho, rhob, neq, nin, ccst, bcst,
     &                    cinf, csup, delta, dwork(pdwo), info)
      ENDIF       
      RETURN
      END
