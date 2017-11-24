c=======================================================================
c
c     subroutine ALLOCVAR                                   
c
c     Asset allocation s.t. maximum Value-at-Risk (upper bound) constraint
c
c     Max[ rho*w]
c      s.t. 
c     VaR(w) <= VaRmax         (Value-at-Risk constraint)
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c
c        w      : portfolio weights 
c        Q      : covariance matrix 
c        rho    : assets performance 
c        VaRmax : VaR upper bound 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocvar ( n, cov, rho,
     &                      neq, nin, ccst, bcst, cinf, csup,
     &                      vrimax, conflp, maxdic, epsdic,
     &                      iwork, dwork,
     &                      wopt, vriopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                            integer
c            cov    : covariance matrix (n*n)                     double
c            rho    : mean returns (n)                            double
c            neq    : number of equality constraints             integer
c            nin    : number of inequality constraints           integer
c            ccst   : constraints matrix (n*(neq + nin))          double
c            bcst   : constraints vector (neq + nin)              double
c            cinf   : lower bound (n)                             double
c            csup   : upper bound (n)                             double
c            vrimax : Value-at-Risk constraint                    double
c            conflp : confidence level parameter ( 0<conflp<1 )   double
c            maxdic : maximum of dichotomy iterations            integer
c            epsdic : precision of dichotomy                      double
c
c     WORKSPACE 
c            iwork  : 3*n + 2*nin + neq + 3                      integer 
c            dwork  : n*(2*n+nin+neq+10) + 2*neq + 4*nin + 4      double
c
c     OUTPUT 
c            wopt   : optimal portfolio (n)                       double
c            vriopt : optimal Value-at-Risk                       double
c            info   : diagnostic argument                        integer
c
c     CALL   
c        YM      : copy a vectorized matrix in a vectorized matrix
c        OPVAR   : computing optimization for ALLOCVAR (cf. OPVAR.F)
c        TESTSDP : test if covariance matrix is SDP
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, maxdic, info
      DOUBLE PRECISION vrimax, conflp, epsdic, vriopt
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(*), 
     &                 bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwo, pdwo, pisdp, pdsdp
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 1.0E-8 )
c
c     external subroutines
      EXTERNAL YM, testsdp, opvar
c     
c     intrinsic functions
      INTRINSIC MAX
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp  : pointer for TESTSDP ( n ) 
      piwo  = 1
c     piwo  : pointer for OPVAR ( 3*n + neq + 2*nin + 3 )
c
c     Total size of iwork array = 3*n + 2*nin + neq + 3
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP ( n*(2*n + 7) )  
      pdwo  = 1
c     pdwo  : pointer for OPVAR
c             ( n*( n + neq + nin + 10 ) + 2*neq + 4*nin + 4 )
c
c     Total size of dwork array = n*(2*n + nin + neq + 10) + 2*neq + 4*nin + 4
c
c-----------------------------------------------------------------------
c
c     test confidence level
      IF ( (conflp .LT. EPS) .or. (conflp .GT. (1. - EPS)) ) THEN
         info = -102
         RETURN
      ENDIF
c
c     covariance matrix SDP test
      CALL testsdp (n , cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF           
c
c     optimization (see OPVAR.F)
      CALL opvar ( n, cov, rho, neq, nin, ccst, bcst, cinf, csup,
     &             vrimax, conflp, maxdic, epsdic,
     &             iwork(piwo), dwork(pdwo),
     &             wopt, vriopt, info )
c
      RETURN
      END
