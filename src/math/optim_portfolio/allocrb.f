c=======================================================================
c
c     subroutine ALLOCRB                            
c
c     Risk Budgeting allocation strategy
c     (Multi-Volatility-constrained)
c
c     Max[ w'*rho - k*w'*Q*w]
c      s.t. 
c     sqrt[w'*cov_i*w] <= sigma(i)  (i-th risk budgeting constraint)
c     C*w <= b                      (general linear constraints) 
c     Cinf <= w <= Csup             (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocrb ( n, cov, Q, kappa, rho,
     &                     neq, nin, ccst, bcst,
     &                     cinf, csup, nbc, class, sigma,
     &                     iwork, dwork, wopt, lambda, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       Q      : matrix (n*n)                                     double
c       kappa  : risk aversion coefficient                        double
c       rho    : expected returns vector (n)                      double
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       ccst   : matrix of constraints (n*(neq+nin))              double
c       bcst   : vector initial of constraints (neq+nin)          double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       nbc    : number of risk budgeting constraints (>=1)      integer 
c       class  : block definition (n)                            integer
c       sigma  : volatilities budgets constraints (nbc)           double 
c
c     WORKSPACE 
c       iwork  : 2*nbc + 5*n + neq + 2*nin + 12                  integer 
c       dwork  : nbc*(nbc+21)/2 + 
c                n*(8*n+20+2*neq+2*nin)+4*neq+6*nin+14            double
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double 
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
      DOUBLE PRECISION kappa
      DOUBLE PRECISION cov(*), Q(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*), sigma(*), lambda(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, npk, totitr, totsim, ncst, infoBFGS, test,
     &        infot, izero
      INTEGER piwbf, pisim, pinfo, pin, pineq, pinin, piclass, piqp,
     &        pydual, pbinf, pbsup, pgrad, pwbf, pdsim, 
     &        pdcov, pdrho, pdcinf, pdcsup, pdsigma, pdc, pdb,
     &        pdcov1, pdQ, pdkappa, pdqp, pdlagr, pdccs, pdbcs
      DOUBLE PRECISION funct, epsbfg, tmp, scal
      DOUBLE PRECISION INFINI, MINFINI, ZERO, EPS, ONE, COEF
      PARAMETER (INFINI=1.E+15,ZERO=0.E0,EPS=1.E-30,ONE=1E0,COEF=-0.5E0)
      PARAMETER (MINFINI = -1.E10)
c
c     external subroutines
      EXTERNAL simrb, simext, testsdp, OVTMCV, IVX, YV, bfgsmvol,
     &         checklinbox, initfeas, IVI
c     
c     intrinsic functions
c
c-----------------------------------------------------------------------
c
c
c     initializations
      info = 0          ! diagnostic argument
      ncst = neq + nin  ! number total of constrains
      izero = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwbf = 1
c     piwbf : pointer for BFGSMVOL who needs (2*nbc + 1)
c     so (2*nbc + 1)
      pisim  = piwbf + (2*nbc + 1)
c     pisim : pointer for simulator SIMRB (5*n + neq + 2*nin + 11)
c     so (2*nbc + 5*n + neq + 2*nin + 12)
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
c     piqp : pointer for SIMRB
c
c     Total size of dwork array = 2*nbc + 5*n + neq + 2*nin + 12 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pydual = 1
c     pydual: pointer for initial value of ydual (nbc)
c     so (nbc)
      pbinf = pydual + ( nbc )
c     pbinf : pointer for inferior bounds binf (nbc)
c     so (2*nbc)
      pbsup = pbinf + ( nbc )
c     pbsup : pointer for superior bounds bsup (nbc)
c     so (3*nbc)
      pgrad = pbsup + ( nbc )
c     pgrad : pointer for gradient at solution (nbc)
c     so (4*nbc)
      pwbf = pgrad + ( nbc )
c     pwbf  : pointer for BFGSMVOL ( nbc*(nbc+11)/2 )   
c     so ( nbc*(nbc+19)/2 )
      pdsim = pwbf + ( nbc*(nbc+11)/2 )
c     pdsim  : pointer for the simulator SIMRB
c             n*(7*n+20+2*neq+2*nin)+4*neq+6*nin+14+nbc
c     so 
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
      pdsigma = pdcsup + ( n )
c     pdsigma : pointer for risk budget constraint (nbc)
      pdc     = pdsigma + ( nbc )
c     pdc : pointer for matrix constraint (n*(neq + nin))      
      pdb     = pdc + ( n*(neq + nin) )
c     pdb : pointer for vector constraint (neq + nin)    
      pdqp    = pdb + ( neq + nin )
c     pdqp : pointer for QP (SIMRB), (n*n + 6*n + 2*nin)
      pdlagr  = pdqp + ( n*n + 6*n + 2*nin )
c     pdlagr : pointer for QP parameter (2*n*n + 3*n + neq + nin)         
      pdcov1  = pdlagr + ( 2*n*n + 3*n + neq + nin )
c     pdcov1  : temporary cov. matrix ( n*n )
      pdccs = pdcov1 + ( n*n )
c     needs n*(neq+nin+1)
      pdbcs = pdccs + n*(neq+nin+1)
c     needs neq+nin+1
c      
c     Total size of dwork array = 4*nbc + ( nbc*(nbc+11)/2 )
c                               + n*(7*n+20+2*neq+2*nin)+4*neq+6*nin+14+nbc (SIMRB)
c                               + n*n + n*(neq+nin+1) + neq+nin+1
c           = nbc*(nbc+21)/2 + n*(8*n+21+3*neq+3*nin)+5*neq+7*nin+15
c
c-----------------------------------------------------------------------
c
c      TODO - test if classes are not all null
c      
c-----------------------------------------------------
c     satisfiability conditions (well-defined problem)
c-----------------------------------------------------
c
c     matrices are SDP ?
      CALL testsdp (n, cov, iwork(pisim), dwork(pdsim), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TESTSDP COV ', info
c      close(unit=1)

      CALL YV ( n*n, Q, dwork(pdQ) )
      IF (kappa .LT. EPS) CALL YV ( n*n, cov, dwork(pdQ) )   
      CALL testsdp (n, dwork(pdQ), iwork(pisim), dwork(pdsim), info)
      IF (info .NE. 0) THEN
         info = -109
         RETURN
      ENDIF

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TESTSDP Q ', info
c      close(unit=1)

c
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neq, nin, ccst, bcst,
     &                iwork(pisim), dwork(pdsim), wopt, info)
      IF (info .EQ. 1111) THEN
        info = -100
        RETURN
      ENDIF

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'INITFEAS 1 ', info
c      close(unit=1)

c
c------------------------------------------------------------
c     construction of work vectors (communication with simrb)
c------------------------------------------------------------
c
c     initialization
      CALL IVX ( nbc, dwork(pydual), ONE )   ! dual
      CALL IVX ( nbc, dwork(pbinf), ZERO )   ! lower bounds
      CALL IVX ( nbc, dwork(pbsup), INFINI ) ! upper bounds
c
c     construct of iwork (communication with simrb)
      iwork(pinfo) = 0       ! diagnostic argument
      iwork(pin)   = n       ! number of asset(s)
      iwork(pineq) = neq     ! number of eq. constraint(s)
      iwork(pinin) = nin     ! number of ineq. constraint(s)
      CALL YVI ( n, class, iwork(piclass) ) ! cov. matrix partition
c
c     construct dwork (communication with simrb)
      CALL YV(n*n, cov, dwork(pdcov))     ! cov. matrix
      CALL YV(n*n, Q, dwork(pdQ))         ! Q matrix
      CALL YV(n, rho, dwork(pdrho))       ! rho     
      dwork(pdkappa) = kappa              ! kappa coef.
      CALL YV(n, cinf, dwork(pdcinf))     ! lower bounds
      CALL YV(n, csup, dwork(pdcsup))     ! upper bounds
      CALL YV(nbc, sigma, dwork(pdsigma)) ! risk budget(s)
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
     &                     + 1.E-4*(MAX(dwork(pdrho + 1),1.E-6))
      ENDIF
c
c     case ncst <> 0 (neq=0, nin=0)
      IF (ncst .NE. 0) THEN
        CALL YM(n, ncst, ccst, dwork(pdc))  ! constraint matrix
        CALL YV(ncst, bcst, dwork(pdb))     ! constraint vector
      ENDIF
c
c-------------------------------------------------------------------
c     Non Linear optimization (BFGS)
c-------------------------------------------------------------------
      epsbfg = 1.E-10   ! precision stop test BFGS
      CALL bfgsmvol ( simrb, simext, nbc, dwork(pydual), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pwbf),
     &                funct, dwork(pgrad), totitr, totsim, info)
     
      CALL YV(nbc, dwork(pydual), lambda) ! ydual -> lambda
      infoBFGS = info ! BFGS diagnostic argument

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'BFGSMVOL ', infoBFGS
c      close(unit=1)

c
c-------------------------------------------------
c     Optimize with optimal dual solution (lambda)
c-------------------------------------------------
c
      CALL IMX ( n, n, dwork(pdcov), ZERO ) ! initialize pdcov
c
c     construct cov. matrix with optimal dual solution
      DO i = 1,nbc
         CALL IVI(n, iwork(piwbf), izero) ! initialize piwbf
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
c     product kappa*Q -> pdcov2     
      CALL PMX(n, n, Q, kappa, dwork(pdcov1))
c
c     kappa*Q + Sum(Cov_i)      
      CALL SM(n, n, dwork(pdcov), dwork(pdcov1), dwork(pdcov))
c
c     linear part: l(i) = -0.5*rho(i)
      scal = -2.
      CALL PVX ( n, rho, COEF, dwork(pdrho) )
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
     &                     + 1.E-4*(MAX(dwork(pdrho + 1),1.E-6))
      ENDIF
c      
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neq, nin, ccst, bcst,
     &                iwork(piqp), dwork(pdqp), wopt, info)
      IF (info .EQ. 1111) THEN
        info = -100
        RETURN
      ENDIF

c     open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'INITFEAS 2', info
c      close(unit=1)

c
c     quadratic solver QP
      CALL qp ( n, dwork(pdcov), dwork(pdrho), neq, nin,
     &          dwork(pdc), dwork(pdb), cinf, csup,
     &          iwork(piqp), dwork(pdqp),
     &          dwork(pdlagr), wopt, info )

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'QUAPRO ', info
c      close(unit=1)

c
c----------------------------------------------------------------------
c     test if the ith-volatility budget constraint is too high/small 
c     (budget constraint not attainable) -> info = 111/110
c----------------------------------------------------------------------
c     
      IF (info .NE. 1001) THEN
        DO i = 1,nbc
            IF (ABS(lambda(i)) .LT. 1.E-15) THEN
                info = 111
            ENDIF
            IF (ABS(lambda(i)) .GT. 1.E+15) THEN
                info = 110
            ENDIF
        ENDDO
      ENDIF
c
c     checking non-emptiness of linear constraints      
c      CALL checklinbox ( n, wopt,
c     &                   neq, nin, ccst, bcst, cinf, csup,
c     &                   dwork(pwbf), infot)
c      IF (infot .LT. 0) THEN
c        info = infot
c        RETURN
c      ENDIF

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'CHECKLINBOX ', info
c      close(unit=1)

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

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TEST 1', info
c      close(unit=1)

         IF (npk .GT. 1) THEN
c
c           initialize temporary cov. matrix 
            CALL IMX ( npk, npk, dwork(pdcov1), ZERO ) 
            CALL IMX ( n, n, dwork(pdcov), ZERO )

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TEST 2', info
c      close(unit=1)

c
c           extracting the i-th block
            CALL YMCPI(n, cov, npk, iwork(piwbf), dwork(pdcov1), infot)

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TEST 3', info
c      close(unit=1)

c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov1), n, iwork(piwbf),
     &                   dwork(pdcov), infot)

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TEST 4', info
c      close(unit=1)

c
c           optimal portfolio volatility: sqrt[wopt'*Cov*wopt]
            CALL OVTMCV(n, dwork(pdcov), wopt, scal)
            scal = SQRT(scal)
            IF (scal .GT. (sigma(i)+1.E-6)) THEN
                info = -120
                RETURN
            ENDIF
         ENDIF
      ENDDO

c      open(unit=1,file='NUMXallocRiskBudgeting.txt',status='unknown')
c      write(1, *) 'TEST 5', info
c      close(unit=1)

c
c     BFGS diagnostic
      IF ((infoBFGS .NE. 0) .AND. (info .EQ. 0)) THEN
        info = infoBFGS
      ENDIF
      RETURN
      END
