c=======================================================================
c
c     subroutine ALLOCITRFR                                     
c
c     Index Tracking allocation srtategy with a risk-free asset
c
c     w(i), i=1,...,n the riskly assets (portfolio weights)
c     w(0), the risk-free asset
c
c     Min[ Tracking Error ]
c      s.t. 
c     (rho*w+rfr*w(0) - index perf. >= delta (relative target return)
c     C*w <= b                               (linear constraints) 
c     Cinf(i) <= w(i) <= Csup(i)             (riskly assets lower/upper bounds)
c     Cinf(0) <= w(0) <= Csup(0)             (risk-free asset lower/upper bounds)
c
c        Q   : covariance matrix
c        rho : expected returns
c        rfr : risk-free rate 
c
c----------------------------------------------------------------------
      SUBROUTINE allocitrfr ( n, cov, rho, covb, rhob, rfr,
     &                        neq, nin, ccst, bcst, cinf, csup, 
     &                        cinfrfr, csuprfr, delta,
     &                        iwork, dwork, wopt, wrf, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : portfolio size                                 integer
c       cov     : covariance matrix (n*n)                         double
c       rho     : expected mean returns (n)                       double
c       covb    : covariance assets-index (n)                     double
c       rhob    : expected mean return of index                   double
c       rfr     : expected risk-free rate                         double
c       neq     : number of equality constraints                 integer
c       nin     : number of inequality constraints               integer
c       ccst    : matrix of constraints (n*(neq+nin))             double
c       bcst    : vector of constraints (neq+nin)                 double
c       cinf    : lower bound (n)                                 double
c       csup    : upper bound (n)                                 double
c       cinfrfr : risk-free lower bound                           double
c       csuprfr : risk-free upper bound                           double
c       delta   : relative target return                          double
c
c     WORKSPACE 
c       iwork   : 3*n + 2*nin + neq + 7                          integer 
c       dwork   : n*(2*neq+2*nin+2*n+23)+6*nin+4*neq+31           double
c                    
c     OUTPUT 
c       wopt    : optimal portfolio (n)                           double
c       wrf     : optimal risk-free asset                         double
c       info    : diagnostic argument                            integer
c
c     CALL   
c       OPIT    : computing optimization for ALLOCIT (cf. OPIT.F)
c       TESTSDP : test if covariance matrix is SDP
c       SEV     : vector sum
c
c-----------------------------------------------------------------------   
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION rhob, delta, rfr, cinfrfr, csuprfr, wrf
      DOUBLE PRECISION cov(*), rho(*), covb(*), cinf(*), csup(*), 
     &                 wopt(*), ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pisdp, piwo, pdsdp, pdwo, pdlagr
      DOUBLE PRECISION sum
c
c     external subroutines
      EXTERNAL testsdp, opitrfr, SEV
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0          ! diagnostic argument
c      
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for workspace of TESTSDP (n)
      piwo  = 1
c     piwo  : pointer for OPITRFR (3*n + 2*nin + neq + 7)
c
c     Total size of iwork array = 3*n + 2*nin + neq + 7 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))  
      pdlagr = 1
c     pdlagr : pointer for lagrange multipliers (n + neq + nin + 3)                
      pdwo  = pdlagr + ( n + neq + nin + 3 )
c     pdwo  : pointer for OPITRFR
c               n*(2*n+2*nin+2*neq+22)+3*neq+5*nin+28
c
c     Total size of dwork array = n*(2*neq+2*nin+2*n+23)+6*nin+4*neq+31
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
c     optimization - traking error minimization (with risk-free asset)
      CALL opitrfr ( n, cov, covb, rho, rhob, rfr, delta,
     &               neq, nin, ccst, bcst, 
     &               cinf, csup, cinfrfr, csuprfr, 
     &               iwork(piwo), dwork(pdwo), wopt, dwork(pdlagr), 
     &               info )
      IF (info .EQ. 1001) THEN
        info = -101
      ENDIF 
c
c     risk-free asset weight
      CALL SEV ( n, wopt, sum )
      wrf = 1.0 - sum
      RETURN
      END
