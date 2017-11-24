c=======================================================================
c
c     subroutine ALLOCSRRFR                                    
c
c     Sharpe Ratio (SR) maximisation strategy with a risk-free asset
c
c     Max[ Sharpe Ratio ]
c      s.t. 
c     C*w <= b                       (linear constraints) 
c     Cinf    <= w    <= Csup        (lower/upper bounds risky assets)
c     Cinfrfr <= wrfr <= Csuprfr     (lower/upper bounds risk-free asset)
c
c     where, Sharpe Ratio := (rho'*w - rfr)/sqrt(w'Qw)
c
c        w    : risky assets weights
c        wrfr : risk-free asset weights 
c        Q    : covariance matrix
c        rho  : assets performance 
c        rfr  : risk-free rate
c
c-----------------------------------------------------------------------
      SUBROUTINE allocsrrfr (n, cov, rho, rfr,
     &                       neq, nin, ccst, bcst, 
     &                       cinf, csup, cinfrfr, csuprfr,
     &                       iwork, dwork, wopt, wrf, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        n       : number of assets                              integer
c        cov     : covariance matrix (n*n)                        double
c        rho     : expected returns (n)                           double
c        rfr     : risk-free rate                                 double
c        neq     : number equality constraints                   integer
c        nin     : number inequality constraints                 integer
c        ccst    : matrix of constraints (n*(neq+nin))            double
c        bcst    : vector initial of constraints (neq+nin)        double
c        cinf    : lower bound (n)                                double
c        csup    : upper bound (n)                                double
c        cinfrfr : risk-free asset lower bound                     double
c        csuprfr : risk-free asset upper bound                     double
c
c     WORKSPACE 
c       iwork    : 7*n + 2*nin + neq + 7                         integer 
c       dwork    : n*(5*n+2*neq+3*nin+35)+4*neq+9*nin+27          double
c
c     OUTPUT 
c       wopt     : optimal portfolio (n)                          double
c       wrf      : optimal risk-free asset                        double
c       info     : diagnostic arguments                          integer
c
c     CALL   
c       TESTSDP : matrix SDP test
c       OPSR    : computing optimization for ALLOCSR
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION rfr, cinfrfr, csuprfr, wrf 
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER p, i, j, neqtot, nintot, ncst, ntot
      INTEGER piwo, pdwo, pdcov, pdrho, pdlagr, pdc, pdb, pdx, 
     &        pdcinf, pdcsup
      DOUBLE PRECISION coef, mu, t, sum, rhomax, ZERO, EPS, NEPS, INF, 
     &                 NINF
      PARAMETER (ZERO = 0.0, EPS = 1.E-12, NEPS = -1.E-12)
      PARAMETER (INF = 1.E+12, NINF = -1.E+12)
c
c     external subroutines
      EXTERNAL testsdp, opmv, qp, YMCPIR, CTSR, IMX, IVX
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0      ! diagnostic argument
      p    = n + 1  ! problem size
      CALL IVX ( n, wopt, ZERO )
      wrf = 0.0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwo  = 1
c     piwo  : pointer for workspaces of TESTSDP, OPMV, QP
c
c     Total size of iwork array = 7*n + 2*nin + neq + 7 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------     
c      pdwo  = 1
c     pdwo  : pointer for workspaces of TESTSDP, OPMV, QP
c             n*( 2*n + neq + nin + 12 ) + 2*neq + 4*nin + 7 
c     pdwo + n*( 2*n + neq + nin + 12 ) + 2*neq + 4*nin + 7 
      pdcov = 1
c     pdcov : pointer for temporary cov. matrix (n+1)*(n+1)
      pdrho = pdcov + ( p*p )
c     pdrho : pointer for temporary rho. verctor (p)
      pdc   = pdrho + ( p )
c     pdc   : pointer for constraint matrix p*(2*n + nin + 2) 
      pdb   = pdc + (  p*(2*n + nin + 2) )
c     pdb   : pointer for constraint vector (2*n + nin + 2)
      pdcinf= pdb + ( 2*n + nin + 2 )
c     pdcinf: pointer for lower bounds (p)
      pdcsup= pdcinf + ( p )
c     pdcsup: pointer for upper bounds (p)
      pdlagr= pdcsup + ( p )
c     pdlagr : pointer for Lagrange multipliers (p + 2 + 2*n + nin)
      pdx   = pdlagr + ( p + 2 + 2*n + nin )
c     pdx    : pointer for quadratic optimal points ( p )
      pdwo  = pdx + p
c     needs 
c           Testsdp needs n*(2*n+7)
c           Quapro needs p*p + 6*p + 2*neq
c           OPMV needs 2*n*(neq+nin+n+9)+4*neq+6*nin+2*n+15
c           so  2*n*(neq+nin+n+9)+4*neq+6*nin+2*n+15

c     Total size of dwork array = n*(2*neq+2*nin+2*n+18)+4*neq+6*nin+2*n+15 + n*(3*n+15+nin) + 12 +3*nin
c                               = n*( 6*n + 2*neq + 3*nin + 35 ) + 4*neq + 9*nin + 27 
c
c-----------------------------------------------------------------------
c
c     case cinfrfr = 100%
      IF ((cinfrfr .GT. 1-EPS).AND.(csuprfr .GT. 1-EPS)) THEN
        wrf = 1.0
        RETURN
      ENDIF
      wrf  = cinfrfr
      coef = 1.0/(1.0 - cinfrfr)
c
c     covariance matrix SDP test
      CALL testsdp (n, cov, iwork(piwo), dwork(pdwo), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     general linear constraints
      ncst   = neq + nin
      neqtot = neq + 1
      ntot   = ncst + 1
      DO i = 1,n
        dwork(pdc + i - 1) =  1.0
      ENDDO
      dwork(pdb) = coef*1.0
c     case ncst = 0 (neq=0, nin=0)
      IF (ncst .NE. 0) THEN
        DO i = 1,n
            DO j = 1,ncst
                dwork(pdc + n + n*(j-1) + i - 1) = ccst(i,j)
            ENDDO
        ENDDO     
        DO j = 1,ncst
            dwork(pdb + j) = coef*bcst(j)
        ENDDO
      ENDIF
c
c     tests constraints compatibility
c      mu = -1.E+8
      DO i = 1,n
        dwork(pdcinf + i - 1) = coef*cinf(i)
        dwork(pdcsup + i - 1) = coef*csup(i)
      ENDDO
c
c     test if max. returns < rfr
      CALL EVMAX ( n, rho, rhomax )
      IF (rhomax .LT. rfr) THEN
        CALL IVX ( n, wopt, ZERO )
        wrf = 1.0
        info = 0
        IF (wrf .GT. (csuprfr+1.E-8)) THEN
            info = -113
        ENDIF
        RETURN
      ENDIF 
c
c     test constraints
      mu = rfr + 1.E-10
      CALL opmv ( n, cov, rho, mu,
     &            neqtot, nin, dwork(pdc), dwork(pdb),
     &            dwork(pdcinf), dwork(pdcsup),
     &            iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .LT. 0) THEN
        info = -100
        RETURN
      ENDIF
c
c     construction of the quadratic part
c     Q = | cov   0 | 
c         | 0 ... 0 |
      CALL IMX ( p, p, dwork(pdcov), ZERO )
      DO i = 1,n
        iwork(piwo + i - 1) = i
      ENDDO
      CALL YMCPIR ( n, cov, p, iwork(piwo), dwork(pdcov), info )
      dwork(pdcov + p*p - 1) = 1.E-12
c
c     initialize linear part to 0
      CALL IVX ( p, dwork(pdrho), ZERO )
c
c     construction of the constraints matrix and vector
      IF ((neq+nin) .GT. 0) THEN
        DO i = 1, (neq+nin)
            dwork(pdwo + i - 1) = coef*bcst(i)
        ENDDO
      ENDIF
      CALL CTSR ( n, rho, rfr, neq, nin, ccst, dwork(pdwo),
     &            dwork(pdcinf), dwork(pdcsup),
     &            neqtot, nintot, dwork(pdc), dwork(pdb))
c
c     initialize lower/upper bounds
      CALL IVX ( p, dwork(pdcinf), ZERO )
      CALL IVX ( p, dwork(pdcsup), INF  )
      dwork(pdcsup + p - 1) = n*INF
c
c     quadratic solver
      CALL qp ( p, dwork(pdcov), dwork(pdrho),
     &          neqtot, nintot, dwork(pdc), dwork(pdb),
     &          dwork(pdcinf), dwork(pdcsup),
     &          iwork(piwo), dwork(pdwo),
     &          dwork(pdlagr), dwork(pdx), info)
      IF (info .NE. 0) THEN
        CALL IVX ( n, wopt, ZERO )
        RETURN
      ENDIF
c
c     re-scale weights
      t = dwork(pdx + p - 1)
c
c     if |t| < EPS then (no solution) 
c     return min. vol. solution            
      IF (ABS(t) .LT. 100*EPS) THEN
c        info = -113
         DO i = 1,n
            wopt(i) = (1.0 - wrf)*wopt(i)
        ENDDO      
        RETURN
      ENDIF
      sum = 0.0
      DO i = 1,n
        wopt(i) = dwork(pdx + i - 1)/t
        sum = sum + wopt(i)
      ENDDO
      DO i = 1,n
        wopt(i) = wopt(i)/sum
      ENDDO
c
c     re-scale result       
      DO i = 1,n
        wopt(i) = (1.0 - wrf)*wopt(i)
      ENDDO
      RETURN
      END
