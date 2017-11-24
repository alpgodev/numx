c=======================================================================
c
c     subroutine ALLOCMVRFR                                    
c
c     Markowitz Mean-Variance strategy with a risk-free asset
c
c     w(i), i=1,...,n the riskly assets (portfolio weights)
c     w(0), the risk-free asset
c
c     Min[ w*Q*w ]
c      s.t. 
c     rho*w + rfr*w(0) >= delta     (target return)
c     C*w <= b                      (linear constraints) 
c     Cinf(i) <= w(i) <= Csup(i)    (riskly assets lower/upper bounds)
c     Cinf(0) <= w(0) <= Csup(0)    (risk-free asset lower/upper bounds)
c
c        Q   : covariance matrix (n by n)
c        rho : expected returns  (n)
c        rfr : risk-free rate    
c
c-----------------------------------------------------------------------
      SUBROUTINE allocmvrfr ( n, cov, rho, rfr,
     &                        neq, nin, ccst, bcst, cinf, csup, 
     &                        cinfrfr, csuprfr, mu,
     &                        iwork, dwork, wopt, wrf, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of assets                               integer
c       cov     : covariance matrix (n*n)                         double
c       rho     : expected returns (n)                            double
c       rfr     : expected risk-free rate                         double
c       neq     : number equality constraints                    integer
c       nin     : number inequality constraints                  integer
c       ccst    : matrix of constraints (nasset*(neq+nin))        double
c       bcst    : vector initial of constraints (neq+nin)         double
c       cinf    : lower bound (n)                                 double
c       csup    : upper bound (n)                                 double
c       cinfrfr : risk-free asset lower bound                     double
c       csuprfr : risk-free asset upper bound                     double
c       mu      : performance target                              double
c
c     WORKSPACE 
c       iwork  : 3*n + 2*nin + neq + 7                           integer 
c       dwork  : n*(2*n+nin+neq+11) + 4*nin + 2*neq + 12          double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                             double
c       wrf    : optimal risk-free asset                           double
c       info   : diagnostic argument                              integer
c
c     CALL   
c       SEV     : vector sum
c       OPMV    : computing optimization for ALLOCMV
c       TESTSDP : matrix SDP test
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION mu, rfr, cinfrfr, csuprfr, wrf
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwo, pisdp, pdwo, pdsdp
      DOUBLE PRECISION sum, mutest
c
c     external subroutines
      EXTERNAL opmvrfr, testsdp, SEV
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0          ! diagnostic argument
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for TESTSDP (n)      
      piwo  = 1
c     piwo  : pointer for OPMVRFR (3*n + 2*nin + neq + 7)
c
c     Total size of iwork array = 3*n + 2*nin + neq + 7 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))              
      pdwo  = 1
c     pdwo  : pointer for OPMVRFR
c             n*(n + neq + nin + 11) + 2*neq  + 4*nin + 12
c
c     Total size of dwork array =  n*(2*n+nin+neq+11) + 4*nin + 2*neq + 12  
c
c-----------------------------------------------------------------------
c
c     covariance matrix SDP (-1.E-50) ?
      CALL testsdp (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     tests constraints compatibility
      mutest = -1.E+8
      CALL  opmvrfr ( n, cov, rho, rfr, mutest, neq, nin, ccst, bcst,
     &                cinf, csup, cinfrfr, csuprfr,
     &                iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -100
        RETURN
      ENDIF
c
c     optimization: mean-variance (with risk-free asset)
      CALL opmvrfr ( n, cov, rho, rfr, mu, neq, nin, ccst, bcst,
     &               cinf, csup, cinfrfr, csuprfr,
     &               iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -101
      ENDIF
c
c     risk-free asset position 
      CALL SEV ( n, wopt, sum )
      wrf = 1.0 - sum
      RETURN
      END
