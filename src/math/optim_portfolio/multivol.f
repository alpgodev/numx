c=======================================================================
c
c     subroutine MULTIVOL                            
c
c     Risk Budgeting allocation strategy 
c
c     Max[ w'*rho ]
c      s.t. 
c     sqrt[w'*Q(i)*w] <= sigma(i)  (volatility constraint)
c     C*w <= b                     (linear constraints) 
c     Cinf <= w <= Csup            (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE multivol ( n, cov, rho, neq, nin, ccst, bcst,
     &                      cinf, csup, nbc, class, sigma,
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
c       nbc    : number of risk budgeting constraints (>=1)      integer 
c       class  : block definition (n)                            integer
c       sigma  : volatilities budgets constraints (nbc)           double 
c
c     WORKSPACE 
c       iwork  : 3*nbc + 4*n + neq + 2*nin + 6                   integer 
c       dwork  : nbc*(nbc+21)/2 + n*(5*n+neq+nin+12)+2*neq+4*nin  double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double 
c       lambda : optimal lambda (nbc)
c       info   : diagnostic argument                             integer
c
c     CALL   
c       SIMMVOL, TESTSDP, OPVOL, OVTMCV, IVX, YV, BFGSMVOL
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simvol, simext
c
c     arguments i/o
      INTEGER n, info, neq, nin, nbc
      INTEGER class(*)
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(n,*), 
     &                 bcst(*), wopt(*), sigma(*), lambda(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, l, k, npk, totitr, totsim, ntot, infoBFGS
      INTEGER piwbf, pisim, pinfo, pin, pineq, pinin, piclass,
     &        pydual, pbinf, pbsup, pgrad, pwbf, psim, 
     &        pdcov, pdrho, pdcinf, pdcsup, pdsigma, pdc, pdb, pdwo,
     &        pdcov1
      DOUBLE PRECISION funct, epsbfg, t, var, volmin, stotal, dualopt
      DOUBLE PRECISION INFINI, MINFINI, ZERO, EPS, ONE, COEF
      PARAMETER (INFINI=1.E15,ZERO=0.E0,EPS=1.E-16,ONE=1E0,COEF=-0.5E0)
      PARAMETER (MINFINI = -1.E10)
c
c     external subroutines
      EXTERNAL testsdp, opvol, OVTMCV, IVX, YV, bfgsmvol
c     
c     intrinsic functions
c
c-----------------------------------------------------------------------
c
c
c     initializations
      info = 0
      ntot = neq + nin
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
c     pydual: pointer for initial value of ydual (nbc)
      pbinf = pydual + ( nbc )
c     pbinf : pointer for inferior bounds binf (nbc)
      pbsup = pbinf + ( nbc )
c     pbsup : pointer for superior bounds bsup (nbc)
      pgrad = pbsup + ( nbc )
c     pgrad : pointer for gradient at solution (nbc)
c
      pwbf = pgrad + ( nbc )
c     pwbf  : pointer for BFGSMVOL ( nbc*(nbc+11)/2 )
c      
      psim = pwbf + ( nbc*(nbc+11)/2 )
c     psim  : pointer for the simulator SIMMVOL
c             ( n*(4*n+neq+nin+12)+nvar+2*neq+4*nin )
      pdcov   = psim
c     pdcov : pointer for cov. matrix (n*n)      
      pdrho   = pdcov + ( n*n )
c     pdrho : pointer for expected returns (n)      
      pdcinf  = pdrho + ( n )
c     pdcinf : pointer for lower bounds (n)      
      pdcsup  = pdcinf + ( n )
c     pdcsup : pointer for upper bounds (n)      
      pdsigma = pdcsup + ( n )
c     pdsigma : pointer for risk budget constraint (nbc)
      pdc     = pdsigma + ( nbc )
c     pdc : pointer for matrix constraint (n*(neq+nin))      
      pdb     = pdc + ( n*(neq + nin) )
c     pdb : pointer for vector constraint (neq+nin)      
      pdwo    = psim + (n + neq + nin + 1)
c     pdwo :       
      pdcov1  = pdb + ( neq + nin )
c     pdcov1  : temporary cov. matrix ( n*n )
c
c     Total size of dwork array = 4*nbc + ( nbc*(nbc+11)/2 )
c                               + n*(4*n+neq+nin+12)+nbc+2*neq+4*nin
c                               + n*n
c           = nbc*(nbc+21)/2 + n*(5*n+neq+nin+12)+2*neq+4*nin
c
c-----------------------------------------------------------------------
c
c      TODO - test if classes are not all null
c
c      open(unit=1,file='MULTIVOL.txt',status='unknown')
c      write(1,*) "--- BEGIN MULTIVOL ---"
c      
c-----------------------------------------------------
c     satisfiability conditions (well-defined problem)
c-----------------------------------------------------
c
c     test if the covariance matrix is SDP
      CALL testsdp (n, cov, iwork(pisim), dwork(psim), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c--------------------------------------------------------------
c     construction of work vectors (communication with simmvol)
c--------------------------------------------------------------
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
      CALL YV ( n, rho, dwork(pdrho) )       ! rho
      CALL YV ( n, cinf, dwork(pdcinf) )     ! lower bounds
      CALL YV ( n, csup, dwork(pdcsup) )     ! upper bounds
      CALL YV ( nbc, sigma, dwork(pdsigma) ) ! risk budget(s)
      CALL YM ( n, ntot, ccst, dwork(pdc) )  ! constraint matrix
      CALL YV ( ntot, bcst, dwork(pdb) )     ! constraint vector
c
c-------------------------------------------------------------------
c     Non Linear optimization (BFGS)
c-------------------------------------------------------------------
      epsbfg = 1.E-10
      CALL bfgsmvol ( SIMMVOL, gesterr, nbc, dwork(pydual), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pwbf),
     &                funct, dwork(pgrad), totitr, totsim, info)
      CALL YV ( nbc, dwork(pydual), lambda ) ! ydual -> lambda
      infoBFGS = info
c      write(1,*) "-------- BFGS info=", info , "-------"
c      write(1,*) "totitr=", totitr ,", totsim= ", totsim 
c      write(1,*) "F(wopt)=", funct
c      write(1,*) "Gradient:"
c      write(1,1) (dwork(pgrad + i - 1), i=1,nbc)
c      write(1,*) "xopt:"
c      write(1,1) (dwork(pydual + i - 1), i=1,nbc)
c      write(1,*) "-----------------------"
c
c     gestion of SIMMVOL errors
c      info = iwork(pinfo)
c      IF (info .LT. 0) RETURN
c
c--------------------------------------------
c     Optimize with optimal lambda
c-------------------------------------------- 
c
      CALL IMX ( n, n, dwork(pdcov), ZERO ) ! initialize pdcov
c
c     construct cov. matrix
      DO i = 1,nbc
c
         CALL IVX ( n, iwork(piwbf), ZERO ) ! initialize piwbf
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

c            write(1,*) "quadratic part, block:",i
c            do l = 1,npk
c              write(1,1)(dwork(pdcov1 + (l-1)*npk + k - 1),k=1,npk)
c            enddo

c
c           product lambda(i)*cov(i)     
            CALL PMX (npk, npk, dwork(pdcov1), lambda(i), dwork(pdcov1))
            
c            write(1,*) "quadratic part * lambda"
c            do l = 1,npk
c              write(1,1)(dwork(pdcov1 + (l-1)*npk + k - 1),k=1,npk)
c            enddo
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov1), n, iwork(piwbf),
     &                   dwork(pdcov), info)
     
c            write(1,*) "Index:"
c            write(1,*) (iwork(piwbf + k - 1), k=1,n)
c            write(1,*) "quadratic part:"
c            do l = 1,n
c              write(1,1)(dwork(pdcov + (l-1)*n + k - 1),k=1,n)
c            enddo

         ENDIF
      ENDDO
c
c     linear part: l(i) = -0.5*rho(i)
      CALL PVX ( n, rho, COEF, dwork(pdrho) )
c
c     quadratic solver
      CALL qp ( n, dwork(pdcov), dwork(pdrho), neq, nin,
     &          ccst, bcst, cinf, csup,
     &          iwork(piclass), dwork(pdb),
     &          dwork(pdcinf), wopt, info )
      IF (info .LT. 0) RETURN
           
c      write(1,*) "QP info=", info
c
c     global risk: SQRT(w'*COV*w)
c      CALL OVTMCV ( n, cov, wopt, var)
c      var = SQRT(var)
c            
c      write(1,*) "wopt=" 
c      write(1,1) (wopt(i), i=1,n)
c      write(1,*) "cov. matrix"
c      do i = 1,n
c         write(1,1)(cov((i-1)*n + j),j=1,n)
c      enddo
c      write(1,*) "cinf="
c      write(1,1) (cinf(i), i=1,n)
c      write(1,*) "csup="
c      write(1,1) (csup(i), i=1,n)
c      write(1,*) "neq=", neq
c      write(1,*) "nin=", nin
c      write(1,*) "constraint matrix:"
c      do i = 1,(neq+nin)
c        write(1,1) (ccst(j,i), j=1,n)
c      enddo
c      write(1,*) "constraint vector:"
c      write(1,1) (bcst(i), i=1,(neq+nin))
c      write(1,*) "linear part:"
c      write(1,1) (dwork(pdrho + i - 1), i=1,n)
c      write(1,*) "quadratic part:"
c      do i = 1,n
c         write(1,1)(dwork(pdcov + (i-1)*n + j - 1),j=1,n)
c      enddo
c    1 format(50f10.6)
c      write(1,*) "--- END MULTIVOL ---"
c      close(unit=1)  

c
c-------------------------------------------------------------   
c     test if the ith-volatility budget constraint is too high/small 
c     (budget constraint not attainable) -> info = 111/110
c-------------------------------------------------------------
c     
      DO i = 1,nbc
        IF (ABS(lambda(i)) .LT. EPS) THEN
            info = 111
        ENDIF
        IF (ABS(lambda(i)) .GT. 1.E+14) THEN
            info = 110
        ENDIF
      ENDDO
c
c     BFGS diagnostic      
      IF ((infoBFGS .NE. 0) .AND. (info .EQ. 0)) THEN
        info = infoBFGS
      ENDIF  
c
c--------------------------------------------
c     Lagrange multiplier (budget constraint)
c-------------------------------------------- 
c
c     TODO    
c
      RETURN
      END
