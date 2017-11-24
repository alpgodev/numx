c=======================================================================
c
c     subroutine MULTIVOLRFR                            
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
      SUBROUTINE multivolrfr ( n, cov, Q, kappa, rho, rfr,
     &                         neq, nin, ccst, bcst,
     &                         cinf, csup, cinfrfr, csuprfr, 
     &                         nbc, class, sigma,
     &                         iwork, dwork, wopt, wrf, lambda, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       Q      : risk aversion matrix (n*n)                       double
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
c              + n*(8*n+2*neq+2*tnin+22)+5*neq+7*tnin+16          double
c              = nbc*(nbc+21)/2 
c              + n*(8*n+2*neq+2*nin+26)+ 5*neq+7*nin +30
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double 
c       wrf    : optimal risk-free asset                          double
c       lambda : optimal lambda (nbc)                             double
c       info   : diagnostic argument                             integer
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
      INTEGER i, j, npk, totitr, totsim, ncst, ntot, tnin, 
     &        infoBFGS, infot, IZERO, test
      INTEGER piwbf, pisim, pinfo, pin, pineq, pinin, piclass, piqp,
     &        pydual, pbinf, pbsup, pgrad, pwbf, pdsim, 
     &        pdcov, pdrho, pdcinf, pdcsup, pdsigma, pdc, pdb,
     &        pdcov1, pdQ, pdkappa, pdcinf0, pdcsup0, pdqp, pdlagr
      DOUBLE PRECISION funct, epsbfg, sum, scal, tmp
      DOUBLE PRECISION INFINI, MINFINI, ZERO, EPS, ONE, COEF, TOL
      PARAMETER (IZERO=0, ZERO=0.E0,EPS=1.E-15,ONE=1.E0,COEF=-0.5E0)
      PARAMETER (TOL=1.E-15, INFINI=1.E+15, MINFINI = -1.E+10)
c
c     external subroutines
      EXTERNAL simmvolrfr, gesterr, testsdp, OVTMCV, IVX, YV,
     &         bfgsmvol, checklinbox, initfeas, IVI, SEV
c     
c     intrinsic functions
c
c-----------------------------------------------------------------------
c
c
c     initializations
      info = 0          ! diagnostic argument
      ntot = neq + nin  ! number total of input constraint
      tnin = nin + 2    ! add risk-free constraints (lower/upper)
      ncst = neq + tnin ! number total of constrains
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwbf = 1
c     piwbf : pointer for BFGSBOX who needs (2*nbc + 1)
      pisim  = piwbf + (2*nbc + 1)
c     pisim : pointer for simulator SIMMVOLRFR (5*n + neq + 2*tnin + 6)
      pinfo = pisim 
c     pinfo : pointer for info (errors of simul), so 1 more
      pin   = pinfo + ( 1 )
c     pin   : pointer for value of n, so 1 more
      pineq = pin + ( 1 )
c     pineq : pointer for value of neq, so 1 more
      pinin = pineq + ( 1 )
c     pinin : pointer for value of nin, so 1 more
      piclass = pinin + ( 1 )
c     piclass: pointer for value of class, so n more
      piqp = piclass + ( n )
c     piqp : pointer for QP (SIMMVOLRFR)  
c
c     Total size of dwork array = (2*nbc + 1) + (5*n + neq + 2*tnin + 6)
c                               = 2*nbc + 5*n + neq + 2*nin + 10
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pydual = 1
c     pydual: pointer for initial value of ydual (nbc)
      pbinf = pydual + ( nbc )
c     pbinf : pointer for inferior bounds binf (nbc)
      pbsup = pbinf + ( nbc )
c     pbsup : pointer for superior bounds bsup (nbc)
      pgrad = pbsup + ( nbc )
c     pgrad : pointer for gradient at solution (nbc)
      pwbf = pgrad + ( nbc )
c     pwbf  : pointer for BFGSMVOL ( nbc*(nbc+11)/2 )
      pdsim = pwbf + ( nbc*(nbc+11)/2 )
c     pdsim  : pointer for the simulator SIMMVOLRFR
c              n*(7*n+2*neq+2*tnin+21)+nbc+4*neq+6*tnin+16
      pdcov   = pdsim
c     pdcov : pointer for cov. matrix (n*n)  
      pdQ     = pdcov + ( n*n )
c     pdQ   : pointer for Q matrix (n*n)      
      pdkappa = pdQ + ( n*n )
c     pdcov : pointer for kappa (1)      
      pdrho   = pdkappa + ( 1 )
c     pdrho : pointer for expected returns (n)      
      pdcinf  = pdrho + ( n )
c     pdcinf : pointer for lower bounds (n)      
      pdcsup  = pdcinf + ( n )
c     pdcsup : pointer for upper bounds (n)
      pdcinf0 = pdcsup + ( n )
c     pdcinf0 : pointer for risk-free lower bounds (1)      
      pdcsup0 = pdcinf0 + ( 1 )
c     pdcsup0 : pointer for risk-free upper bounds (1)      
      pdsigma = pdcsup0 + ( 1 )
c     pdsigma : pointer for risk budget constraint (nbc)
      pdc     = pdsigma + ( nbc )
c     pdc : pointer for matrix constraint (n*(neq + tnin))      
      pdb     = pdc + ( n*(neq + tnin) )
c     pdb : pointer for vector constraint (neq + tnin)    
      pdqp    = pdb + ( neq + tnin )
c     pdqp : pointer for QP (SIMMVOLRFR), n*(neq+tnin+2*n+15)+2*neq +4*tnin+13
      pdlagr  = pdqp + ( n*(5*n+neq+tnin+18)+3*neq+5*tnin+13 )
c     pdlagr : pointer for QP Lagrange parameter (n + neq + tnin)           
      pdcov1  = pdlagr + ( n + neq + tnin )
c     + n*n
c
c     Total size of dwork array = 4*nbc + ( nbc*(nbc+11)/2 )
c                               + n*(7*n+2*neq+2*tnin+21)+nbc+4*neq+6*tnin+16 (SIMMVOLRFR)
c                               + n*n
c                               + n+neq+tnin
c              = nbc*(nbc+21)/2 + n*(8*n+2*neq+2*tnin+22)+5*neq+7*tnin+16
c-----------------------------------------------------------------------
c
c      TODO - test if classes are not all null
c      
c-----------------------------------------------------
c     satisfiability conditions (well-defined problem)
c-----------------------------------------------------
c
c     test if the covariance matrix is SDP
      CALL testsdp (n, cov, iwork(pisim), dwork(pdsim), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF    
      IF (kappa .LT. EPS) THEN 
        CALL YV ( n*n, cov, dwork(pdQ) )
      ELSE     
        CALL YV ( n*n, Q, dwork(pdQ) )
        CALL testsdp (n, dwork(pdQ), iwork(pisim), dwork(pdsim), info)
        IF (info .NE. 0) THEN
            info = -109
            RETURN
        ENDIF
      ENDIF  
c      
c     case sigma(1)=0 
      IF (nbc .EQ. 1) THEN
        IF (ABS(sigma(1)) .LT. EPS) THEN
            CALL IVX ( n, wopt, ZERO )
            IF (csuprfr .GT. (1.0-10*EPS)) THEN
                wrf = 1.0
            ELSE
                wrf  = cinfrfr
                info = 110
            ENDIF
            RETURN
        ENDIF      
      ENDIF  
c
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neq, nin, ccst, bcst,
     &               iwork(piqp), dwork(pdqp), wopt, info )
      IF (info .EQ. 1111) THEN
        info = -100
        RETURN
      ENDIF  
c
c-----------------------------------------------------------------
c     construction of work vectors (communication with simmvolrfr)
c-----------------------------------------------------------------
c
c     initialization
      CALL IVX ( nbc, dwork(pydual), ONE )   ! dual
      CALL IVX ( nbc, dwork(pbinf), ZERO )   ! lower bounds
      CALL IVX ( nbc, dwork(pbsup), INFINI ) ! upper bounds
c
c     construct of iwork (communication with simmvol)
      iwork(pinfo) = 0       ! diagnostic argument
      iwork(pin)   = n       ! number of asset(s)
      iwork(pineq) = neq     ! number of eq. constraint(s)
      iwork(pinin) = nin     ! number of ineq. constraint(s)
      CALL YVI ( n, class, iwork(piclass) ) ! cov. matrix partition
c
c     construct dwork (communication with simmvol)
      CALL YV ( n*n, cov, dwork(pdcov) )     ! cov. matrix
      CALL YV ( n*n, Q, dwork(pdQ) )         ! Q matrix (risk aversion)
      CALL YV ( n, rho, dwork(pdrho) )       ! rho
      DO i = 1,n
        dwork(pdrho + i - 1) = dwork(pdrho + i - 1) - rfr
      ENDDO
      dwork(pdkappa) = kappa                 ! kappa coef.
      CALL YV ( n, cinf, dwork(pdcinf) )     ! lower bounds
      CALL YV ( n, csup, dwork(pdcsup) )     ! upper bounds
      CALL YV ( 1, cinfrfr, dwork(pdcinf0) ) ! risk-free lower bound
      CALL YV ( 1, csuprfr, dwork(pdcsup0) ) ! risk-free upper bound
      CALL YV ( nbc, sigma, dwork(pdsigma) ) ! risk budget(s)
c
c     case rho(1) = rho(2) = ... = rho(n)
      tmp = rho(1)
      test = 0
      DO i = 2,n
       IF (ABS(tmp - rho(i)) .GT. EPS) THEN
        test = 1
       ENDIF
      ENDDO
      IF (test .EQ. 0) THEN
        dwork(pdrho + 1) = dwork(pdrho + 1) 
     &                   + 1.E-4*(MAX(dwork(pdrho + 1),1.E-6))
      ENDIF
c
c     case ncst <> 0 (neq=0, nin=0)
      IF (ntot .NE. 0) THEN
        CALL YM ( n, ntot, ccst, dwork(pdc) )  ! constraint matrix
        CALL YV ( ntot, bcst, dwork(pdb) )    ! constraint vector
      ENDIF
c
c     linear constraints (C*w <= b)
      DO i = 1,n
        dwork(pdc + ntot*n + i - 1)     =  1.0 
        dwork(pdc + (ntot+1)*n + i - 1) = -1.0
      ENDDO
      dwork(pdb + ntot)     = 1.0 - cinfrfr
      dwork(pdb + ntot + 1) = csuprfr - 1.0
c
c-------------------------------------------------------------------
c     Non Linear optimization (BFGS)
c-------------------------------------------------------------------
      epsbfg = 1.E-10   ! precision stop test
      CALL bfgsmvol ( simmvolrfr, gesterr, nbc, dwork(pydual), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pwbf),
     &                funct, dwork(pgrad), totitr, totsim, info)
      CALL YV ( nbc, dwork(pydual), lambda ) ! ydual -> lambda
      infoBFGS = info ! BFGS diagnostic argument
c
c-------------------------------------------------
c     Optimize with optimal dual solution (lambda)
c-------------------------------------------------
c
      CALL IMX ( n, n, dwork(pdcov), ZERO ) ! initialize pdcov
c
c     construct cov. matrix with optimal dual solution
      DO i = 1,nbc
         CALL IVI(n, iwork(piwbf), IZERO) ! initialize piwbf
c
c        dimension of block i
         npk = 0
         DO j = 1,n
            IF (class(j) .EQ. i) THEN
               npk = npk + 1
               iwork(piwbf + npk - 1) = j
            ENDIF
         ENDDO
         IF (npk .GT. 1) THEN
c
c           initialize temporary cov. matrix 
            CALL IMX ( npk, npk, dwork(pdcov1), ZERO ) 
c         
c           extract block i
            CALL YMCPI ( n, cov, npk, iwork(piwbf), dwork(pdcov1), info)
c
c           product lambda(i)*cov(i)     
            CALL PMX (npk, npk, dwork(pdcov1), lambda(i), dwork(pdcov1))
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov1), n, iwork(piwbf),
     &                   dwork(pdcov), info)
         ENDIF
      ENDDO
c            
c     product kappa*Q -> pdcov1     
      CALL PMX (n, n, Q, kappa, dwork(pdcov1))
c
c     kappa*Q + Sum(Cov_i)
      CALL SM ( n, n, dwork(pdcov), dwork(pdcov1), dwork(pdcov))      
c
c     linear part: l(i) = -0.5*rho(i)
      CALL PVX ( n, rho, COEF, dwork(pdrho) )
c
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neq, tnin, dwork(pdc), dwork(pdb),
     &                iwork(piqp), dwork(pdqp), wopt, info)
      IF (info .EQ. 1111) THEN
        info = -100
        RETURN
      ENDIF  
c
c     quadratic solver 
      CALL qp ( n,dwork(pdcov), dwork(pdrho), neq, tnin,
     &          dwork(pdc), dwork(pdb), cinf, csup,
     &          iwork(piqp), dwork(pdqp),
     &          dwork(pdlagr), wopt, info )
c
c     risk-free asset
      CALL SEV ( n, wopt, sum )
      wrf = 1.0 - sum
c
c-------------------------------------------------------------   
c     test if the ith-volatility budget constraint is too high/small 
c     (budget constraint not attainable) -> info = 111/110
c-------------------------------------------------------------
c     
      IF (info .NE. 1001) THEN
        DO i = 1,nbc
            IF (ABS(lambda(i)) .LT. 1.E-15) THEN
                info = 111
            ENDIF
            IF (ABS(lambda(i)) .GT. 1.E+14) THEN
                info = 110
            ENDIF
        ENDDO
      ENDIF
c
c     checking non-emptiness of linear constraints      
      CALL checklinbox ( n, wopt, neq, tnin,
     &                   dwork(pdc), dwork(pdb),
     &                   cinf, csup,
     &                   dwork(pwbf), infot )
      IF (infot .LT. 0) THEN
        info = infot
        RETURN
      ENDIF
c
c     construct cov. matrix with optimal dual solution
      DO i = 1,nbc
         DO j = 1,n
            iwork(piwbf + i - 1) = 0
         ENDDO
c
c        dimension of the i-th block
         npk = 0
         DO j = 1,n
            IF (class(j) .EQ. i) THEN
               npk = npk + 1
               iwork(piwbf + npk - 1) = j
            ENDIF
         ENDDO
         IF (npk .GT. 1) THEN
c
c           initialize temporary cov. matrix 
            CALL IMX ( npk, npk, dwork(pdcov1), ZERO ) 
            CALL IMX ( n, n, dwork(pdcov), ZERO ) 
c         
c           extracting the i-th block
            CALL YMCPI(n, cov, npk, iwork(piwbf), dwork(pdcov1), infot)
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov1), n, iwork(piwbf),
     &                   dwork(pdcov), infot)        
c
c           variance of the i-th block, w'*Cov*w     
            CALL OVTMCV(n, dwork(pdcov), wopt, scal)
            scal = SQRT(scal)
            IF (scal .GT. (sigma(i)+1.E-4)) THEN
                info = -120
                RETURN
            ENDIF
         ENDIF
      ENDDO   
c
c     BFGS diagnostic (very complex case)     
      IF ((infoBFGS .NE. 0) .AND. (info .EQ. 0)) THEN
        info = infoBFGS
      ENDIF  
      RETURN
      END
