c=======================================================================
c
c     subroutine ALLOCMV                                     
c
c     Markowitz Mean-Variance strategy
c
c     Min[ w*Q*w ]
c      s.t. 
c     rho*w >= delta           (target return)
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocmv ( n, cov, rho,
     &                     neq, nin, ccst, bcst, cinf, csup, mu,
     &                     iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of assets                               integer
c       cov     : covariance matrix (n*n)                         double
c       rho     : expected returns (n)                            double
c       neq     : number equality constraints                    integer
c       nin     : number inequality constraints                  integer
c       ccst    : matrix of constraints (n*(neq+nin))             double
c       bcst    : vector initial of constraints (neq+nin)         double
c       cinf    : lower bound (n)                                 double
c       csup    : upper bound (n)                                 double
c       mu      : performance target                              double
c
c     WORKSPACE 
c       iwork   : 3*n+2*nin+neq+6                                integer 
c       dwork   : 2*n*(neq+nin+n+9)+4*neq+6*nin+2*n+15            double
c
c     OUTPUT 
c       wopt    : optimal portfolio (n)                           double
c       info    : diagnostic argument                            integer
c
c     CALL   
c       OPMV    : computing optimization for ALLOCMV
c       TESTSDP : test if covariance matrix is SDP
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION mu
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwo, pdwo, pisdp, pdsdp
c
c     external subroutines
      EXTERNAL opmv, testsdp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for TESTSDP (n)      
      piwo  = 1
c     piwo  : pointer for OPMV 3*n+2*nin+neq+6 
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 6 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))              
      pdwo  = 1
c     pdwo  : pointer for OPMV
c            2*n*(neq+nin+n+9)+4*neq+6*nin+2*n+15 
c
c     Total size of dwork array = n*(2*n + 2*nin + 2*neq + 18) + 4*neq + 
c                                                 6*nin + 15
c
c-----------------------------------------------------------------------
c
c     covariance matrix SDP ?
      CALL testsdp (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     optimization - Mean/Variance (cf. OPMV.F)
      CALL opmv ( n, cov, rho, mu,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -101
      ENDIF 
      RETURN
      END
