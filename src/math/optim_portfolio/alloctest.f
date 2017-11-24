c=======================================================================
c
c     subroutine ALLOCTEST
c
c     Max. parameters for Allocation strategies 
c       - Max. target return
c       - Max. relative targer return
c       - Min. volatility
c       - Min. VaR
c       - Min. CVaR
c
c-----------------------------------------------------------------------
      SUBROUTINE alloctest ( n, cov, rho, rhob, alpha,
     &                       neq, nin, ccst, bcst, cinf, csup,  
     &                       iwork, dwork, 
     &                       targetmax, deltamax,
     &                       volmin, varmin, cvarmin,
     &                       info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n         : portfolio size                               integer
c       cov       : covariance matrix (n*n)                       double
c       rho       : mean returns vector (n)                       double
c       rhob      : mean return of index                          double
c       alpha     : confidence level parameter (0<alpha<1)        double
c       neq       : number equality constraints                  integer
c       nin       : number inequality constraints                integer
c       ccst      : matrix of constraints (nasset*(neq+nin))      double
c       bcst      : vector initial of constraints (neq+nin)       double
c       cinf      : lower bound (n)                               double
c       csup      : upper bound (n)                               double
c
c     WORKSPACE 
c       iwork     : 4*n + 2*nin + neq + 1                        integer 
c       dwork     : n*(4*n + 16) + 3*nin + neq                    double
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
      DOUBLE PRECISION rhob, alpha, targetmax, deltamax, volmin, varmin,
     &                 cvarmin
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(*),bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, piw, pdcov, pdrho, pdw, pdlagr, pdwopt, infotmp,
     &        pisdp, pdsdp
      DOUBLE PRECISION a, sum, mean, p, q, zc, ZERO, PI, EPS
      PARAMETER (ZERO = 0.D0, PI = 3.1415926535898, EPS = 1.D-15)
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL IMX, IVX, YV, qp, OVTMCV, EXARET, testsdp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw  = 1
c     piw  : pointer for QP workspaces (3*n + 2*nin + neq + 1)
      pisdp = piw + ( 3*n + 2*nin + neq + 1 )
c     pisdp  : pointer for TESTSDP (n)  
c
c     Total size of iwork array = 4*n + 2*nin + neq + 1
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcov  = 1
c     pdcov  : pointer for cov. matrix(n*n)
      pdrho  = pdcov + ( n*n )
c     pdrho  : pointer for vector rho, so (n) more
      pdw    = pdrho +( n )
c     pdw    : pointer for QP workspaces (n*n + 6*n + 2*nin)
      pdlagr = pdw + ( n*n + 6*n + 2*nin )
c     pdlagr : pointer for Lagrange multiplier, (n + nin + neq)
      pdwopt = pdlagr + ( n + nin + neq )
c     pdwopt : pointer for wopt, (n)
      pdsdp = pdwopt + ( n )
c     pdsdp : pointer for TESTSDP (n*(2*n + 7)) 
c
c     Total size of dwork array = n*n 
c                               + n
c                               + n*n + 6*n + 2*nin
c                               + n + nin + neq
c                               + n
c                               + n*(2*n + 7)
c                              =  n*(4*n + 16) + 3*nin + neq
c
c-----------------------------------------------------------------------
c
c     initialize quadratic part to zero
      CALL IMX ( n, n, dwork(pdcov), ZERO )
c
c     initialize linear part to -rho
      a = -1.0
      CALL PVX ( n, rho, a, dwork(pdrho) )
c
c     quadratic solver: min[-rho'*w]
      CALL qp ( n, dwork(pdcov), dwork(pdrho),
     &          neq, nin, ccst, bcst, cinf, csup,
     &          iwork(piw), dwork(pdw),
     &          dwork(pdlagr), dwork(pdwopt), infotmp )
      IF (infotmp .LT. 0) info = infotmp
c
c     max. target
      sum = 0.0
      DO i = 1,n
        sum = sum + rho(i)*dwork(pdwopt+i-1)
      ENDDO
      targetmax = sum
      deltamax = targetmax - rhob
c      
c     initialize linear part to -rho      
      CALL IVX ( n, dwork(pdrho), ZERO )
c
c     covariance matrix SDP test
      CALL testsdp (n, cov, iwork(pisdp), dwork(pdsdp), infotmp)
      IF (infotmp .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c      
c     quadratic solver: min[w'*cov*w]    
      CALL qp ( n, cov, dwork(pdrho),
     &          neq, nin, ccst, bcst, cinf, csup,
     &          iwork(piw), dwork(pdw),
     &          dwork(pdlagr), dwork(pdwopt), infotmp )
      IF (infotmp .LT. 0) info = infotmp
c
c     sqrt[w'*cov*w] (portfolio volatility)
      CALL OVTMCV ( n, cov, dwork(pdwopt), volmin)
      volmin = SQRT(volmin)
c
c     mean return
      CALL EXARET ( n, rho, dwork(pdwopt), mean, infotmp )
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
c     computing Normal -VaR := max[0, min(1, -zc*volatility - mean)] 
      varmin = - volmin * zc + mean
      varmin = MIN(varmin, 1.)
      varmin = MAX(varmin, -1.)
c
c     coeff := [sqrt(2*pi)*exp[(zc/2)^2]*(1-alpha)]^(-1)
      zc = SQRT(2.0*PI)*EXP(((zc*zc)/4.0))*(1.0 - alpha)
      zc = 1.0/zc
c     
c     computing Normal Conditional Value-at-Risk
c     CVaR := max[0, min(1, -coeff*volatility - mean)] 
      cvarmin = - volmin * zc + mean
      cvarmin = MIN(cvarmin, 1.)
      cvarmin = MAX(cvarmin, -1.)
      RETURN
      END
