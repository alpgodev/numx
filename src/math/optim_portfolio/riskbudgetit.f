c=======================================================================
c
c     subroutine RISKBUDGETIT Index Tracking                            
c
c     Risk Budgeting allocation strategy (Index Tracking constraints)
c
c     Max[ w'*rho -kappa(1/2*w'*Gamma*w-Covb*w]
c      s.t. 
c     sqrt[(w-w_b)'*Gamma_i*(w-w_b)] <= sigma(i)**2-var_bench+w_b'*Gamma_i*w_b  (volatility constraint)
c     C*w <= b                     (linear constraints) 
c     Cinf <= w <= Csup            (lower/upper bounds)
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE riskbudgetit ( n, cov, kappa, rho, covb, varb, neq,
     &                          nin, ccst, bcst, cinf, csup,
     &                          nbc, class, sigma,
     &                          iwork, dwork,
     &                          wopt, lambda, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : expected returns vector (n)                      double
c       covb   : covariance assets-index (n)                      double
c       varb   : variance of the benchmark                        double
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
c       iwork  : 3*nbc + 17*n + neq + 2*nin + 12                 integer 
c       dwork  : nbc*(nbc+25)/2+n*(13*n+53+2*neq
c                   +2*nin)+6*neq+8*nin+14                        double
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
      INTEGER n, info, neq, nin, nbc
      INTEGER class(*)
      DOUBLE PRECISION cov(*), rho(*), covb(*), varb, cinf(*), csup(*)
      DOUBLE PRECISION ccst(n,*), bcst(*), wopt(*), sigma(*), lambda(*)
      DOUBLE PRECISION  kappa
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     Local Variables
c       
      double precision kappa2, scal, scal2, dzero, epsil
      
      integer pdQ, pdrhoIT, pdbcstIT, pdcinfIT, pdcsupIT, pdsigmaIT
      integer pdomegab, pdw, pdwalloc, pdwoptIT, pdcov, pdcov1
      
      integer piw, piwalloc
      
      integer iosort, infot, i, j, npk
      iosort = 6
      
      pdQ = 1
c     needs n*n
c      
      
      pdomegab = pdQ + n*n
c     needs n
c       
      pdrhoIT = pdomegab + n
c     needs n
c      so n*n+2*n
      pdbcstIT = pdrhoIT + n
c     needs neq + nineq
c      so n*n+2*n + neq + nineq
      pdcinfIT = pdbcstIT + neq + nin
c     needs n 
c      so n*n+3*n + neq + nineq
      pdcsupIT = pdcinfIT + n
c     needs n 
c      so n*n+4*n + neq + nineq   
      pdsigmaIT = pdcsupIT + n
c     needs nbc
c      so n*n+4*n + neq + nineq + nbc
      pdwoptIT = pdsigmaIT + nbc
c     needs n 
c      so n*n+5*n + neq + nineq + nbc
      pdw = pdwoptIT + n
c     needs n*(4*n+28) + nbc + nineq + neq
c      so    n*n+n*(5*n+33) + 2*(nbc + nineq + neq)
      pdwalloc = pdw + n*(4*n+28) + nbc + nin + neq
c     needs  nbc*(nbc+21)/2 + 
c           n*n+n*(8*n+20+2*neq+2*nin)+4*neq+6*nin+14
c
c     so nbc*(nbc+25)/2+n*(14*n+53+2*neq+2*nin)+6*neq+8*nin+14
      pdcov = pdwalloc 
c     needs n*n
c      
      pdcov1 = pdcov + n*n
c     needs n*n
c
      
      piw = 1
c     needs 12*n + nbc 
c
      piwalloc = piw + 12*n + nbc  
c     needs 2*nbc + 5*n + neq + 2*nin + 12  
c          so 3*nbc + 17*n + neq + 2*nin + 12
      call buildrbit(n, cov, kappa, rho, covb, varb, 
     &   neq, nin, ccst, bcst, cinf, csup, nbc, class, sigma,  
     &   iwork(piw), dwork(pdw), dwork(pdomegab), dwork(pdrhoIT), 
     &   dwork(pdbcstIT), dwork(pdcinfIT), dwork(pdcsupIT), 
     &   dwork(pdsigmaIT), info)
      if (info .LT. 0) then
        return
      endif 
      kappa2 = kappa/2

      
      call ym(n, n, cov, dwork(pdQ))
      call allocrb ( n, cov, dwork(pdQ), kappa2, dwork(pdrhoIT),
     &               neq, nin,
     &               ccst, dwork(pdbcstIT),
     &               dwork(pdcinfIT), dwork(pdcsupIT),
     &               nbc, class, dwork(pdsigmaIT),
     &               iwork(piwalloc), dwork(pdwalloc),
     &               dwork(pdwoptIT), lambda, info)
      call SV(n, dwork(pdwoptIT), dwork(pdomegab), wopt)
      if (info .eq. 1001) then
        info = -101
        return
      endif 
c
c     building cov. matrix with optimal dual solution
      DO i = 1,nbc
c       
         DO j = 1,n
            iwork(piwalloc + i - 1) = 0
         ENDDO
c         CALL IVX ( n, iwork(piwbf), ZERO ) ! initialize piwbf
c
c        dimension of block i
         epsil = 1.e-8
         dzero = 0.
         npk = 0
         DO j = 1,n
            IF (class(j) .EQ. i) THEN
               npk = npk + 1
               iwork(piwalloc + npk - 1) = j
            ENDIF
         ENDDO
         IF (npk .GT. 1) THEN
c
c           initialize temporary cov. matrix 
            CALL IMX ( npk, npk, dwork(pdcov1), dzero) 
            CALL IMX ( n, n, dwork(pdcov), dzero) 
c         
c           extract block i
            CALL YMCPI(n, cov, npk, iwork(piwalloc), dwork(pdcov1)
     &     , infot)
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov1), n, iwork(piwalloc),
     &                   dwork(pdcov), infot)
c
c           optimal portfolio variance: w'*Cov*w          
            CALL OVTMCV(n, cov,  dwork(pdwoptIT), scal)
            IF (scal .GT. (dwork(pdsigmaIT+i-1)**2+epsil)) THEN
                info = -120
                RETURN
            ENDIF
            CALL OVTMCV(n, cov,  wopt, scal)
            CALL XV(n, wopt, covb, scal2) 
            IF (scal .gt. (sigma(i)**2+2*scal2-varb+epsil)) THEN
                info = -120
                RETURN
            ENDIF
         ENDIF
      ENDDO 
      RETURN
      END
