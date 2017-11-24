c=======================================================================
c
c     subroutine ALLOCMV1                                     
c
c     This function implements the Markowitz strategy (Mean-Variance)
c
c     Min[ w*Q*w ]
c      s.t. 
c     rho*w >= delta           (target return)
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c     #assets < pmax           (max. nb. assets) 
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE ALLOCMV1 ( n, cov, rho,
     &                     neq, nin, ccst, bcst, cinf, csup, mu, pmax,
     &                     iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                             integer
c            cov    : covariance matrix (n*n)                     double
c            rho    : mean returns vector (n)                     double
c            neq    : number equality constraints                integer
c            nin    : number inequality constraints              integer
c            ccst   : matrix of constraints (nasset*(neq+nin))    double
c            bcst   : vector initial of constraints (neq+nin)     double
c            cinf   : lower bound (n)                             double
c            csup   : upper bound (n)                             double
c            mu     : performance target                          double
c            pmax   : max. number of assets                      integer 
c
c     WORKSPACE 
c            iwork  : 3*n + 2*nin + neq + 3                       integer 
c            dwork  : n*(2*n+nin+neq+9) + 4*nin + 2*neq + 4        double
c
c     OUTPUT 
c            wopt   : optimal portfolio (n)                        double
c            info   : diagnostic argument                         integer
c
c     CALL   
c       YM      : copy a vectorized matrix in a vectorized matrix
c       OPMV    : computing optimization for ALLOCMV (cf. OPMV.F)
c       TESTSDP : test if covariance matrix is SDP
c       EVMAX   :        
c       EVMIN   :
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin, pmax
      DOUBLE PRECISION mu
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, piwo, pdwo, pisdp, pdsdp, nassets, ind
      DOUBLE PRECISION mutest, wmin, EPS
      PARAMETER (EPS = 1.E-8)
c
c     external subroutines
      EXTERNAL OPMV, TESTSDP
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
c     piwo  : pointer for OPMV (3*n + 2*nin + neq + 3)
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))              
      pdwo  = 1
c     pdwo  : pointer for OPMV
c            ( n*( n + neq + nin + 9 ) + 2*neq + 4*nin + 4  )
c
c     Total size of dwork array = n*(2*n + nin + neq + 9) + 4*nin + 2*neq + 4
c
c-----------------------------------------------------------------------
c
c     test if pmax < n
      IF (pmax .GT. n) THEN
         info = -111
         RETURN
      ENDIF      
c
c     covariance matrix SDP ?
      CALL TESTSDP (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     tests constraints compatibility
      mutest = -1.E+8
      CALL OPMV ( n, cov, rho, mutest,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info )           
      IF (info .EQ. 1001) THEN
        info = -100
        RETURN
      ENDIF  
      GOTO 999
c
c     optimization (cf. OPMV.F)
999   CALL OPMV ( n, cov, rho, mu,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -101
        RETURN
      ENDIF
c
c     nb. assets lower/upper constraints
      nassets = n
      DO i = 1,n
        IF (wopt(i) .LT. EPS) THEN
            wopt(i) = 101
            cinf(i) = 0
            csup(i) = 0
            nassets = nassets - 1
        ENDIF
      ENDDO
c
c     if nassets > pmax then goto 999
      IF (nassets .GT. pmax) THEN
        CALL EVMINI ( n, wopt, wmin, ind )
        cinf(ind) = 0
        csup(ind) = 0
        GOTO 999
      ENDIF
c
c     update small wopt(i)      
      DO i =1,n
        IF (wopt(i) .GT. 100) THEN
            wopt(i) = 0           
        ENDIF      
      ENDDO       
c
      RETURN
      END
