c=======================================================================
c
c     subroutine MULTIVOLX (Expert Version)                            
c
c     This function implements multi-Volatility-constrained optimization
c
c     Max[ w'*rho ]
c      s.t. 
c     w'*Q*w <= sigma2         (volatility constraint)
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE multivolx (n,cov,rho, neq, nin, ccst, bcst, cinf, csup,
     &                      nbc, class, sigma,
     &                      dxmin, df1, epsabs, mode, hessian, 
     &                      totitr, totsim,
     &                      iwork, dwork, wopt, lambda, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : expected returns vector (n)                      double
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       ccst   : matrix of constraints (nasset*(neq+nin))         double
c       bcst   : vector initial of constraints (neq+nin)          double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       nbc    : number of risk budgeting constraints (>1)       integer 
c       class  : block definition (n)                            integer
c       sigma  : volatilities budgets constraints (nbc)           double 
c       dxmin  : precision on variables stop criter (nbc)         double
c       df1    : criter decrease at the first iteration           double
c       epsabs : precision stop test                              double
c       mode   : initialization mode                             integer 
c       hessian: initial Hessian (nbc*nbc)                        double             
c       totitr : number max. of iterations                       integer
c       totsim : number max. of simulator call                   integer
c
c     WORKSPACE 
c       iwork  : 3*nbc + 4*n + neq + 2*nin + 6                   integer 
c       dwork  : nbc*(nbc+19)/2 + n*(5*n+neq+nin+12)+2*neq+4*nin  double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double 
c       lambda : optimal lambda (nbc)
c       info   : diagnostic argument                             integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simmvol, gesterr
c
c     arguments i/o
      INTEGER n, info, neq, nin, nbc, mode, totitr, totsim
      INTEGER class(*)
      DOUBLE PRECISION epsabs, df1
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(n,*), 
     &                 bcst(*), wopt(*), sigma(*), lambda(*), dxmin(*),
     &                 hessian(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, npk, ntot, infoBFGS
      INTEGER piwbf, pisim, pinfo, pin, pineq, pinin, piclass,
     &        pydual, pbinf, pbsup, pgrad, pwbf, psim, 
     &        pdcov, pdrho, pdcinf, pdcsup, pdsigma, pdc, pdb, pdwo,
     &        pdcov1
      DOUBLE PRECISION funct, t, var, volmin, stotal
      DOUBLE PRECISION INFINI, MINFINI, ZERO, EPS, ONE, COEF
      PARAMETER (INFINI=1.E10,ZERO=0.E0,EPS=1.E-50,ONE=1E0,COEF=-0.5E0)
      PARAMETER (MINFINI = -1.E10) 
c
c     external subroutines
      EXTERNAL testsdp, opvol, OVTMCV, IVX, YV, bfgsmvolx, qp
c     
c     intrinsic functions
c
c-----------------------------------------------------------------------
c
c
c     initializations
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwbf = 1
c     piwbf : pointer for BFGSBOX who needs (2*nbc + 1)
      pisim  = piwbf + (2*nbc + 1)
c     pisim : pointer for simulator SIMMVOL (5*n + neq + 2*nin + 5)
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
c
c     Total size of dwork array = (2*nbc + 1) + (5*n + neq + 2*nin + 5)
c                               = 2*nbc + 5*n + neq + 2*nin + 6 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pydual = 1
c     pydual: pointer for initial value of ydual, so nbc more
      pbinf = pydual + ( nbc )
c     pbinf : pointer for inferior bounds binf, so nbc more
      pbsup = pbinf + ( nbc )
c     pbsup : pointer for superior bounds bsup, so nbc more
      pgrad = pbsup + ( nbc )
c     pgrad : pointer for gradient at solution, so nbc more
c
      pwbf = pgrad + ( nbc )
c     pwbf  : pointer for BFGSMVOL ( nbc*(nbc+9)/2 )
c      
      psim = pwbf + ( nbc*(nbc+9)/2 )
c     psim  : pointer for the simulator SIMMVOL
c             ( n*(4*n+neq+nin+12)+nvar+2*neq+4*nin )
      pdcov   = psim
      pdrho   = pdcov + ( n*n )
      pdcinf  = pdrho + ( n )
      pdcsup  = pdcinf + ( n )
      pdsigma = pdcsup + ( n )
      pdc     = pdsigma + ( nbc )
      pdb     = pdc + ( n*(neq + nin) )
      
      pdwo    = psim + (n + neq + nin + 1)
      
      pdcov1  = pdb + ( neq + nin )
c     pdcov1  : temporary cov. matrix ( n*n )
c
c     Total size of dwork array = 4*nbc + ( nbc*(nbc+9)/2 )
c                               + n*(4*n+neq+nin+12)+nbc+2*neq+4*nin
c                               + n*n
c           = nbc*(nbc+19)/2 + n*(5*n+neq+nin+12)+2*neq+4*nin
c
c-----------------------------------------------------------------------
c
      IF (nbc .LE. 0) RETURN
c
c      open(unit=1,file='MULTIVOL.txt',status='unknown')
c      write(1,*) "--- BEGIN MULTIVOL ---"
c      
c-----------------------------------------------------
c     satisfiability conditions (well-defined problem)
c-----------------------------------------------------
c
c     covariance matrix SDP test
      CALL testsdp (n, cov, iwork(pisim), dwork(psim), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     test if global budget volatility >= min. volatility
      t = -1.E+15
c     optimization (cf. OPVOL.F)
      CALL opvol ( n, cov, rho, t,
     &             neq, nin, ccst, bcst, cinf, csup,
     &             iwork(pisim), dwork(pdwo),
     &             wopt, dwork(psim), info)
     
c      write(1,*) "OPVOL info=", info
c      do i = 1,n
c        write (1,*) wopt(i)
c      enddo
c
c     if (w'*Q*w) < volb2 -> global volatility budget const. too small !
      CALL OVTMCV ( n, cov, wopt, var)
      volmin = SQRT(var)
      stotal = sigma(1)
      
c      write(1,*) "----TEST -109----"
c      write(1,*) "Min. Volatility=", volmin
c      write(1,*) "Global Risk Budget=", sigma(1)
c      write(1,*) "wopt:" 
c      do i = 1,n
c        write(1,*) wopt(i)
c      enddo
c      write(1,*) "---------------"
      
c      write(1,*) "min. volatility=", volmin
c      write(1,*) "global risk budget=", smin 
      
      IF (stotal .LT. volmin) THEN
        info = -109
        RETURN
      ENDIF 
c--------------------------------------------------------------
c     construction of work vectors (communication with simmvol)
c--------------------------------------------------------------
c     initialization of ydual
      CALL IVX ( nbc, dwork(pydual), ONE )
c
c     initializations of bounds
      CALL IVX ( nbc, dwork(pbinf), ZERO )
      CALL IVX ( nbc, dwork(pbsup), INFINI )
c
c     construct of iwork (communication with simmvol)
      iwork(pinfo) = 0
      iwork(pin)   = n
      iwork(pineq) = neq
      iwork(pinin) = nin
      CALL YVI ( n, class, iwork(piclass) )
c
c     construct dwork (communication with simmvol)      
c     cov matrix: cov
      CALL YV ( n*n, cov, dwork(pdcov) )
c     expected return: rho
      CALL YV ( n, rho, dwork(pdrho) )
c     lower bounds
      CALL YV ( n, cinf, dwork(pdcinf) )
c     upper bounds
      CALL YV ( n, csup, dwork(pdcsup) )
c     risk budgets
      CALL YV ( nbc, sigma, dwork(pdsigma) )
c     constraint matrix
      ntot = neq + nin
      CALL YM ( n, ntot, ccst, dwork(pdc) )
c     constraint vector
      CALL YV ( ntot, bcst, dwork(pdb) )
c     
c-------------------------------------------------------------------
c     optimization BFGS (Expert Version)
c
c     (the second parameter is not used here, we can replace GESTERR
c      by an other name of subroutine)
c-------------------------------------------------------------------
      CALL bfgsmvolx ( simmvol, gesterr, nbc, dwork(pydual),
     &                 dxmin, df1, epsabs, mode, totitr, totsim, 
     &                 dwork(pbinf), dwork(pbsup), 
     &                 iwork(piwbf), dwork(pwbf),
     &                 funct, dwork(pgrad), info )
c
c     stock BFGS info
      infoBFGS = info
c
c      DEBUG TRACE
c      write(1,*) "-------- BFGS info=", info , "-------"
c      write(1,*) "totitr=", totitr ,", totsim= ", totsim 
c      write(1,*) "F(wopt)=", funct
c      write(1,*) "Gradient:"
c      write(1,1) (dwork(pgrad + i - 1), i=1,nbc)
c      write(1,*) "xopt:"
c      write(1,1) (dwork(pydual + i - 1), i=1,nbc)
c      write(1,*) "-----------------------"
      
      IF ( info .LT. 0 ) THEN
c        close(unit=1)
        RETURN
      ENDIF
c
c     gestion of SIMMVOL errors
c      info = iwork(pinfo)
c      IF (info .LT. 0) RETURN
c
c--------------------------------------------
c     Optimize with optimal lambda
c-------------------------------------------- 
c
c     copy ydual -> lambda
      CALL YV ( nbc, dwork(pydual), lambda )
c     
c 
      CALL IMX ( n, n, dwork(pdcov1), ZERO )
c
c     construct cov. matrix
      IF (nbc .GT. 1) THEN
      DO i = 2,nbc
c
c        dimension of block i 
         npk = 0
         DO j = 1,n
            IF (class(j) .EQ. i-1) THEN
               npk = npk + 1
               iwork(piwbf + npk - 1) = j
            ENDIF
         ENDDO
         IF (npk .GT. 1) THEN
c
c           initialize temporary cov. matrix 
            CALL IMX ( npk, npk, dwork(pdc), ZERO ) 
c         
c           extract block i
            CALL YMCPI ( n, cov, npk, iwork(piwbf),
     &                   dwork(pdc), info )
c
c           product lambda(i)*cov(i)     
            CALL PMX (npk, npk, dwork(pdc), lambda(i), dwork(pdc))
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdc), n, iwork(piwbf),
     &                   dwork(pdcov1), info)
c           YMCPIR : info = 0 allways by construction
         ENDIF
      ENDDO
      ENDIF
c
c     global risk budget
      CALL PMX (n, n, cov, lambda(1), dwork(pdcov))   
c
c     sum      
      CALL SM( n, n, dwork(pdcov), dwork(pdcov1), dwork(pdcov))
c
c     linear part: l(i) = -0.5*rho(i)
      CALL PVX ( n, rho, COEF, dwork(pdrho) )
c
c     quadratic solver
      CALL qp ( n,dwork(pdcov), dwork(pdrho), ccst, bcst,
     &          cinf, csup, neq, nin,
     &          iwork(piclass), dwork(pdb),
     &          dwork(pdcinf), wopt, info )
           
c      write(1,*) "QP info=", info
c
c     global risk: SQRT(w'*COV*w)
      CALL OVTMCV ( n, cov, wopt, var)
      var = SQRT(var)
      
c      write(1,*) "Global risk=", var
      
c      write(1,*) (class(i), i = 1,n)
      
c     block risk constraints
      IF (nbc .GT. 1) THEN
      DO i = 2,nbc
      
c         write(1,*) "Class", i
c
c        dimension of block i 
         npk = 0
         DO j = 1,n
            IF (class(j) .EQ. (i-1)) THEN
               npk = npk + 1
               iwork(piwbf + npk - 1) = j
            ENDIF
         ENDDO
         
c         write(1,*) "Nb Class", npk
         
         IF (npk .GT. 1) THEN
c
c           initialize cov1 
            CALL IMX ( npk, npk, dwork(pdc), ZERO )   
            CALL IMX ( n, n, dwork(pdc + npk*npk), ZERO )   
c         
c           extract block i
            CALL YMCPI ( n, cov, npk, iwork(piwbf),
     &                   dwork(pdc), info )
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdc), n, iwork(piwbf),
     &                   dwork(pdc + npk*npk), info)
c
c           SQRT[w'*Q(i)*w]          
            CALL OVTMCV ( n, dwork(pdc + npk*npk), wopt, var)
            IF (var .LT. EPS) THEN
                var = 0.0
            ELSE
                var = SQRT(var)
            ENDIF      
c            write(1,*) "risk(",i,")=", var
         ENDIF
      ENDDO
      ENDIF
      
c      write(1,*) "wopt:" 
c      do i = 1,n
c        write(1,*) wopt(i)
c      enddo
c      write(1,*) "cinf="
c      do i = 1,n
c        write(1,1) cinf(i)
c      enddo
c      write(1,*) "csup="
c      do i = 1,n
c        write(1,1) csup(i)
c      enddo
c      write(1,*) "neq=", neq
c      write(1,*) "nin=", nin
c      write(1,*) "constraint matrix:"
c      do i = 1,(neq+nin)
c        write(1,1) (ccst(j,i), j=1,n)
c      enddo
c      write(1,*) "constraint vector:"
c      do i = 1,(neq+nin)
c        write(1,1) bcst(i)
c      enddo
c      write(1,*) "linear part:"
c      do i = 1,n
c        write(1,1) dwork(pdrho + i - 1)
c      enddo
c      write(1,*) "quadratic part:"
c      do i = 1,n
c         write(1,1)(dwork(pdcov + (i-1)*n + j - 1),j=1,n)
c      enddo
c    1 format(50f10.6)
c      write(1,*) "--- END MULTIVOL ---"
c      close(unit=1)  

c
c---------------------------------------------------------   
c     test if the volatility budget constraint is too high 
c     (budget constraint not attainable) -> info = 112
c---------------------------------------------------------
c     
      IF (lambda(1) .LT. EPS) THEN
        info = 112
      ENDIF  
      
      IF (infoBFGS .NE. 0) THEN
        info = infoBFGS
      ENDIF  
c
c--------------------------------------------
c     Lagrange multiplier (budget constraint)
c-------------------------------------------- 
c
      RETURN
      END
