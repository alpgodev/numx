c=======================================================================
c
c     subroutine ALLOCTESTRFR
c
c     Max. parameters for Allocation strategies with risk-free asset
c       - Max. target return
c       - Max. relative targer return
c       - Min. volatility
c       - Min. VaR
c       - Min. CVaR
c
c-----------------------------------------------------------------------
      SUBROUTINE alloctestrfr ( n, cov, rho, rfr, rhob, alpha,
     &                          neq, nin, ccst, bcst, 
     &                          cinf, csup, cinfrfr, csuprfr,  
     &                          iwork, dwork, 
     &                          targetmax, deltamax, 
     &                          volmin, varmin, cvarmin,
     &                          info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n         : portfolio size                               integer
c       cov       : covariance matrix (n*n)                       double
c       rho       : mean returns vector (n)                       double
c       rfr       : risk-free rate                                double
c       rhob      : mean return of index                          double
c       alpha     : confidence level parameter (0<alpha<1)        double
c       neq       : number equality constraints                  integer
c       nin       : number inequality constraints                integer
c       ccst      : matrix of constraints (nasset*(neq+nin))      double
c       bcst      : vector initial of constraints (neq+nin)       double
c       cinf      : lower bound (n)                               double
c       csup      : upper bound (n)                               double
c       cinfrfr   : risk-free asset lower bound                   double
c       csuprfr   : risk-free asset upper bound                   double
c
c     WORKSPACE 
c       iwork     : 6*n + 2*nin + neq + 13                       integer 
c       dwork     : n*(9*n+2*neq+2*nin+27)+7*nin+5*neq+41         double
c
c     OUTPUT 
c       targetmax : maximum target return (MV)                    double
c       deltamax  : maximum relative target return (IT)           double
c       volmin    : minimum volatility                            double
c       varmin    : minimum Value-at-Risk (VaR)                   double
c       cvarmin   : minimum Conditional Value-at-Risk (CVaR)      double
c       info      : diagnostic argument                          integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION rfr, rhob, alpha, targetmax, deltamax, volmin, 
     &                 varmin, cvarmin, cinfrfr, csuprfr
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(*),bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, nbc, infotmp
      INTEGER piw, piclass, pdw, pdwopt, pdQ
      DOUBLE PRECISION mu, sigma, lagr, wrf, kappa, p, q, zc,
     &                 targetmin, EPS, ZERO, PI, DELTA
      PARAMETER (ZERO = 0.D0, PI = 3.1415926535898)
      PARAMETER (EPS = 1.E-15, DELTA = 1.E-8)
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL allocmvrfr, allocrbrfr, OVTMCV, IMX, EXARET
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
      targetmax = 0.0 
      deltamax  = 0.0 
      volmin    = 0.0
      varmin    = 0.0
      cvarmin   = 0.0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piclass = 1
c     piclass : pointer for asset claasification (n)      
      piw     = piclass + ( n )
c     piw     : workspace 2*nbc + 5*n + 2*nin + neq + 11
c     Total size of iwork array = 6*n + 2*nin + neq + 13
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwopt = 1
c     pdwopt : pointer for wopt ( n )
      pdQ    = pdwopt + ( n )
c     pdQ    : pointer for cov. matrix (n*n)    
      pdw    = pdQ + ( n*n )
c     pdw    : workspace 
c                nbc*(nbc+21)/2 
c              + n*(8*n+2*neq+2*nin+26)+ 7*nin+5*neq+30
c              = n*(8*n+2*neq+2*nin+26)+ 7*nin+5*neq+41
c     Total size of dwork array = n*(9*n+2*neq+2*nin+27)+7*nin+5*neq+41  
c
c-----------------------------------------------------------------------     
c
c     min. volatility     
      mu = -1.0E+15         
      CALL allocmvrfr ( n, cov, rho, rfr,
     &                  neq, nin, ccst, bcst,
     &                  cinf, csup, cinfrfr, csuprfr, mu,
     &                  iwork(piw), dwork(pdw),
     &                  dwork(pdwopt), wrf, infotmp )
      IF (infotmp .LT. 0) info = infotmp
c
c     sqrt[w'*cov*w] (portfolio volatility)
      CALL OVTMCV ( n, cov, dwork(pdwopt), volmin)
      IF (volmin .GT. EPS) THEN
        volmin = SQRT(volmin) + DELTA
      ELSE
        volmin  = 0.0
      ENDIF
      CALL EXARET(n,rho,dwork(pdwopt),mu,infotmp)
      IF (infotmp.LT.0) RETURN
      targetmin = mu + wrf*rfr - DELTA
c
c     max. return: max[rho'*w]
      kappa = 0.0
      CALL IMX ( n, n, dwork(pdQ), ZERO )
      nbc = 1
      DO i = 1,n
        iwork(piclass + i - 1) = 1
      ENDDO
      sigma =  1.0E+15
      CALL allocrbrfr ( n, cov, dwork(pdQ), kappa, rho, rfr,
     &                  neq, nin, ccst, bcst,
     &                  cinf, csup, cinfrfr, csuprfr,
     &                  nbc, iwork(piclass), sigma,
     &                  iwork(piw), dwork(pdw),
     &                  dwork(pdwopt), wrf, lagr,
     &                  infotmp )
      IF (infotmp .LT. 0) info = infotmp
c      
c     min return
      CALL EXARET ( n, rho, dwork(pdwopt), mu, infotmp )
      targetmax = mu + wrf*rfr 
      deltamax = targetmax - rhob
c
c     test if: 0 < alpha < 1
      IF ((alpha .LT. EPS).OR.(alpha .GT. (1.0-EPS))) THEN
        info = -102
        RETURN
      ENDIF
c      
c     computing inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p = alpha
      q = 1. - p
      zc = dinvnr(p,q)
c
c     case zc=0
      IF (ABS(zc) .GT. EPS) THEN
        zc = 1.0/zc
      ELSE
        varmin  = 0.0
        cvarmin = 0.0
        RETURN
      ENDIF
c     
c     computing Normal -VaR := max[0, min(1, -zc*volatility - mean)] 
      varmin = - volmin * zc + targetmin
      varmin = MIN(varmin, 1.)
      varmin = MAX(varmin, -1.)
c
c     coeff := [sqrt(2*pi)*exp[(zc/2)^2]*(1-alpha)]^(-1)
      zc = SQRT(2.0*PI)*EXP(((zc*zc)/4.0))*(1.0 - alpha)
c     
c     computing Normal Conditional Value-at-Risk
c     CVaR := max[0, min(1, -coeff*volatility - mean)] 
      cvarmin = - volmin * zc + targetmin
      cvarmin = MIN(cvarmin, 1.)
      cvarmin = MAX(cvarmin, -1.)
      RETURN
      END
