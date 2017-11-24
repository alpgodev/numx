c=======================================================================
c
c     subroutine ALLOCRBRFR 
c
c     Risk Budgeting allocation strategy with a risk-free asset
c
c     w(i), i=1,...,n the riskly assets (portfolio weights)
c     w(0), the risk-free asset
c
c     Max[ w'*rho - k*w'*Q*w]
c      s.t. 
c     sqrt[w'*cov_i*w] <= sigma(i)  (i-th risk budgeting constraint)
c     C*w <= b                      (general linear constraints) 
c     Cinf(i) <= w(i) <= Csup(i)    (riskly assets lower/upper bounds)
c     Cinf(0) <= w(0) <= Csup(0)    (risk-free asset lower/upper bounds)
c
c        Q   : covariance matrix
c        rho : expected returns
c        rfr : risk-free rate 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocrbrfr ( n, cov, Q, kappa, rho, rfr,
     &                        neq, nin, ccst, bcst,
     &                        cinf, csup, cinfrfr, csuprfr, 
     &                        nbc, class, sigma,
     &                        iwork, dwork, wopt, wrf, lambda, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       Q      : matrix (n*n)                                     double
c       kappa  : risk aversion coefficient                        double
c       rho    : expected returns vector (n)                      double
c       rfr    : expected risk-free rate                          double
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       ccst   : matrix of constraints (nasset*(neq+nin))         double
c       bcst   : vector initial of constraints (neq+nin)          double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       cinfrfr: risk-free lower bound                            double
c       csuprfr: risk-free upper bound                            double
c       nbc    : number of risk budgeting constraints (>=1)      integer 
c       class  : block definition (n)                            integer
c       sigma  : volatilities budgets constraints (nbc)           double 
c
c     WORKSPACE 
c       iwork  : 2*nbc + 5*n + neq + 2*tnin + 7                  integer 
c              = 2*nbc + 5*n + neq + 2*nin + 11
c       dwork  : nbc*(nbc+21)/2 
c              + n*(8*n+2*neq+2*nin+26)+5*neq+7*nin+30            double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double 
c       wrf    : optimal risk-free asset                          double
c       lambda : optimal lambda (nbc)                             double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin, nbc
      INTEGER class(*)
      DOUBLE PRECISION wrf, rfr, kappa, cinfrfr, csuprfr
      DOUBLE PRECISION cov(*), Q(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*), sigma(*), lambda(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, infotmp, iter, piw, pdw
      DOUBLE PRECISION ncinfrfr, ncsuprfr, coef, sum, vol, EPS, ZERO
      PARAMETER (ZERO=0.E0, EPS=1.E-15)
c
c     external subroutines
      EXTERNAL multivolrfr, PVX2, SEV, testlinearcst
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0          ! diagnostic argument
      iter = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1    
c
c     Total size of dwork array = 2*nbc + 5*n + neq + 2*nin + 11
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdw = 1
c
c     Total size of dwork array = nbc*(nbc+21)/2 
c                               + n*(8*n+3*neq+3*nin+28) 
c                               + 5*neq + 7*nin + 31
c
c-----------------------------------------------------------------------
c    
c     case sigma(1)=0 
      IF (nbc .EQ. 1) THEN
        IF (ABS(sigma(1)) .LT. EPS) THEN
            CALL IVX ( n, wopt, ZERO )
            IF (csuprfr .GT. (1.0-1.E-8)) THEN
                wrf = 1.0
            ELSE
                wrf  = cinfrfr
                info = 110
            ENDIF
            RETURN
        ENDIF      
      ENDIF  
c
      CALL multivolrfr ( n, cov, Q, kappa, rho, rfr,
     &                   neq, nin, ccst, bcst,
     &                   cinf, csup, cinfrfr, csuprfr, 
     &                   nbc, class, sigma,
     &                   iwork(piw), dwork(pdw), 
     &                   wopt, wrf, lambda, info )  
      infotmp = info
      GOTO 10
c
c     case info = -120 (risk budget cst. not repected)
 10   iter = iter + 1
      IF (iter .GT. 2) GOTO 200
      IF (info .EQ. -120) THEN
        IF (nbc .EQ. 1) THEN
c           portfolio volatility
            info = 0    
            CALL OVTMCV(n, cov, wopt, vol)
            IF (vol .GT. 1.E-15) THEN
                vol = SQRT(vol)
                coef = sigma(1)/vol
            ELSE 
                coef = 0.0
            ENDIF         
            CALL PVX2 ( n, wopt, coef )
            CALL SEV ( n, wopt, sum )
            wrf = 1.0 - sum
            IF ((wrf .LT. (cinfrfr - 1.E-8)) .OR. 
     &          (wrf .GT. (csuprfr + 1.E-8))) THEN
                info = 110
                GOTO 200
            ENDIF 
            CALL testlinearcst ( n, wopt, neq, nin, ccst, bcst, cinf,
     &                           csup, dwork(pdw), info )
            IF (info .LT. 0) THEN
                info = -115
                GOTO 50
            ENDIF
            GOTO 200    
        ELSE
            GOTO 60
        ENDIF
       
 50     coef = 0.0 ! risk-free rate = 0%
        CALL multivolrfr ( n, cov, Q, kappa, rho, coef,
     &                     neq, nin, ccst, bcst,
     &                     cinf, csup, cinfrfr, csuprfr,
     &                     nbc, class, sigma,
     &                     iwork(piw), dwork(pdw),
     &                     wopt, wrf, lambda, info )
        GOTO 10
       
 60     IF (cinfrfr .EQ. 0) THEN
            ncinfrfr = 0.0
            ncsuprfr = 0.0
            GOTO 100
        ENDIF
        IF (csuprfr .LT. 1) THEN
            ncinfrfr = csuprfr
            ncsuprfr = csuprfr
            GOTO 100
        ENDIF
100     CALL multivolrfr ( n, cov, Q, kappa, rho, rfr,
     &                     neq, nin, ccst, bcst,
     &                     cinf, csup, ncinfrfr, ncsuprfr,
     &                     nbc, class, sigma,
     &                     iwork(piw), dwork(pdw),
     &                     wopt, wrf, lambda, info )
        IF (info .EQ. -120) THEN
            info = 110
        ENDIF 
      ENDIF  
      GOTO 200
200   IF ((info.NE.1001).AND.(info.NE.-108).AND.(info.NE.-109)) THEN
        DO i = 1,nbc
            IF (ABS(lambda(i)) .LT. 1.E-14) THEN
                info = 111
            ENDIF
            IF (ABS(lambda(i)) .GT. 1.E+14) THEN
                info = 110
            ENDIF
        ENDDO
      ENDIF
      IF ((info .EQ. 0).AND.(infotmp.NE.-120).AND.(infotmp.NE.0)) THEN
        info = infotmp
      ENDIF
      
      RETURN
      END
