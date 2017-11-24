c=======================================================================
c
c     subroutine ALLOCVOLRFR                            
c
c     This function implements Volatility-constrained optimization with
c     risk-free rate
c
c     Max[ w'*rho ]
c      s.t. 
c     w'*Q*w <= sigma2         (volatility constraint)
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c     Cinf <= w(0) <= Csup 
c
c        w   : portfolio weights
c        Q   : covariance matrix
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocvolrfr ( n, cov, rho,
     &                         neq, nin, ccst, bcst, cinf, csup, 
     &                         volb, volp, rfr, cinfrfr, csuprfr, 
     &                         iwork, dwork, wopt, lambda, info )
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
c       volb   : volatility budget (>0)                           double 
c       volp   : precision of the volatility cst. (>0)            double
c       rfr    : risk-free rate                                   double 
c       cinfrfr: risk-free rate lower bound                       double
c       csuprfr: risk-free rate upper bound                       double
c
c     WORKSPACE 
c       iwork  : 3*n + 2*nin + neq + 7                           integer 
c       dwork  : n*( 2*n+neq+nin+13 ) + 4*nin + 2*neq + 12        double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n+1)                          double
c       lambda : lagrange multiplier corresponding to the risk 
c                budget constraint (shadow cost)                  double 
c       info   : diagnostic argument                             integer
c
c     CALL   
c       EVMAX   : max. value of a vector
c       EVMIN   : min. value of a vector 
c       OVTMCV  : quadratic form: x'*A*x -> scalar
c       YM      : copy a vectorized matrix in a vectorized matrix
c       OPVOL   : computing optimization for ALLOCVOL (cf. modules/OPVOL.F)
c       TESTSDP : SDP test of a matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION epskat, volb, volp, lambda, rfr, cinfrfr, csuprfr
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, niter, maxiter,
     &        pisdp, piwo, pdsdp, pdwo, pdlagr, pdrho, pdwopt
      DOUBLE PRECISION eps, rhomax, rhomin, rhomean, tmin, tmax, t, sum,
     &                 error, mincinf, var, volmin, volopt,epsilon
      PARAMETER (eps = 1.E-50)
c
c     external subroutines
      EXTERNAL EVMAX, EVMIN, YM, opvolrfr, testsdp, OVTMCV
c     
c     intrinsic functions
      INTRINSIC MAX, MIN, ABS
c
c-----------------------------------------------------------------------
c
c     initializations 
      info  = 0      ! diagnostic argument
      maxiter = 1000 ! max. number of iterations
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for TESTSDP (n)      
      piwo  = 1
c     piwo  : pointer for internal workspaces of OPVOL
c             needs( 3*n + 2*nin + neq + 7 )
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 7 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))  
      pdrho = 1
c     pdrho : pointer for excess return ( n )
      pdwopt = pdrho + ( n )
c     pdwopt : pointer for wopt vector ( n )                        
      pdwo  = pdwopt + ( n )
c     pdwo  : pointer for internal workspaces of OPVOLRFR
c             needs( n*( n + neq + nin + 10 ) 
c                    + neq + 3*nin + 9  )
      pdlagr = pdwo + ( n*(n+neq+nin+10) + neq + 3*nin + 9 )
c     pdlagr : pointer for Lagr of OPVOLRFR
c              needs ( n + neq + nin + 3 )      
c
c     Total size of dwork array = 2*n
c                               + n*( n + neq + nin + 10 )
c                               + neq + 3*nin + 9
c                               + n + neq + nin + 3
c     = n*( 2*n + neq + nin + 13 ) + 4*nin + 2*neq + 12
c
c
c-----------------------------------------------------------------------
c      open(unit=1,file='ALLOCVOLRFR.txt',status='unknown')
c      write(1,*) "Status=", 'OK'
c      do i = 1,n     
c        do j = 1,ncsttot
c            write(1,*) ccstot(i,j)
c        enddo
c      enddo
c       write(1,*) "-----------------------"
c      close(unit=1) 
c      RETURN
c
c-----------------------------------------------------
c     satisfiability conditions (well-defined problem)
c-----------------------------------------------------
c
c     covariance matrix SDP test
      CALL testsdp (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c      
c     excess expected return: rho(i)-rfr
      DO i = 1,n
        dwork(pdrho + i - 1) = rho(i) - rfr
      ENDDO
c
c     test if budget volatility >= min. volatility
      t = -1.E+15
c     optimization (cf. rne/modules/OPVOLRFR.F)
      CALL opvolrfr ( n, cov, dwork(pdrho), t,
     &                neq, nin, ccst, bcst, cinf, csup,
     &                cinfrfr, csuprfr,
     &                iwork(piwo), dwork(pdwo),
     &                dwork(pdwopt), dwork(pdlagr), info)
      sum = 0.0
      DO i = 1,n
        wopt(i+1) = dwork(pdwopt+i-1)
        sum = sum + wopt(i+1)
      ENDDO
      wopt(1) = 1. - sum 
c
c     if (w'*Q*w) < volb2 -> volatility budget const. too small !
      CALL OVTMCV ( n, cov, dwork(pdwopt), var)
      volmin = SQRT(var)
      IF (volb .LT. volmin) THEN
        info = -109
        RETURN
      ENDIF 
c
c     the risk-budget constraint must be saturated. 
c     a sufficient condition is that all the expected returns are not null 
c     (case rho(i) = 0)
      DO i = 1,n
        IF (ABS(rho(i)) .LT. eps) THEN
            info = -110
            RETURN
        ENDIF
      ENDDO
      
c------------------------------------
c     initialization of tmin/tmax
c------------------------------------
c
c     min. lower bound
      CALL EVMIN(n, cinf, mincinf)
c
c     case long positions only
      IF (mincinf .GE. 0) THEN
c
c       max./min. return(s)
        CALL EVMAX(n, rho, rhomax)
        rhomax = MAX(ABS(rfr),rhomax) 
        CALL EVMIN(n, rho, rhomin)
        rhomin = MIN(ABS(rfr),rhomin)
        tmax = rhomax
        tmin = rhomin
c
c     case long/short positions        
      ELSE  
        rhomax = ABS(rfr)
        DO i = 1,n
            rhomax = MAX(rhomax, ABS(rho(i))) 
        ENDDO
        tmax = ABS(rfr)*(MAX(ABS(cinfrfr), ABS(csuprfr)))
        DO i = 1,n
            tmax = tmax + rhomax * (MAX(ABS(cinf(i)), ABS(csup(i))))  
        ENDDO
        tmin = -tmax
      ENDIF
c
c     myzero := max(1.E-12, mean(rho))
      CALL MV ( n, rho, rhomean )
      epsilon = MAX(1.E-12, rhomean*1.E-8)
c
c      open(unit=1,file='ALLOCVOL.txt',status='unknown')
c
c     dichotomy loop
      error = 1.E+15
      niter = 0
      IF ( ABS(tmax - tmin) .LT. epsilon ) THEN
        info = 112
        RETURN
      ENDIF 
c      
c---------------------------------
c     dichotomy on tagret perf.
c---------------------------------
c
      DO WHILE ((niter .LT. maxiter) .AND. (error .GT. volp))
c 
        niter = niter + 1
        t = ((tmin + tmax)/2.0) - rfr    
c
c       optimization with risk-free rate (cf. OPVOLRFR.F)     
        CALL opvolrfr ( n, cov, dwork(pdrho), t,
     &                  neq, nin, ccst, bcst, cinf, csup,
     &                  cinfrfr, csuprfr,
     &                  iwork(piwo), dwork(pdwo),
     &                  dwork(pdwopt), dwork(pdlagr), info )
c
c       sqrt[w'*Q*w] (portfolio volatility)
        CALL OVTMCV ( n, cov, dwork(pdwopt), var)      
        volopt = SQRT(var)
c
c       if sqrt(w'*Q*w) > volb, then tmax = t, else tmin = t       
        IF ((volopt .GT. volb) .OR. (info .NE. 0)) THEN
            tmax = t
        ELSE
            tmin = t
        ENDIF 
c             
c       if |tmax - tmin| < eps, then EXIT
        IF ( ABS(tmax - tmin) .LT. epsilon ) THEN
            info = 112
            RETURN
        ENDIF        
c
c       error = |sqrt(w'*Q*w) - volb|     
        error = ABS(volopt - volb)
c        
c      open(unit=1,file='ALLOCVOL.txt',status='unknown')
c       write(1,*) "mu=", t
c       write(1,*) "Iteration=", niter
c      do i = 1,n     
c        do j = 1,ncsttot
c            write(1,*) ccstot(i,j)
c        enddo
c      enddo
c       write(1,*) "-----------------------"
c      close(unit=1) 
c      
      ENDDO
      
      IF (niter .EQ. maxiter) THEN
        info = 111
      ENDIF
c
c--------------------------------------------
c     optimal portfolio: w=[w(0),...,w(n)]
c--------------------------------------------
      sum = 0.0
      DO i = 1,n
        wopt(i+1) = dwork(pdwopt + i - 1)
        sum = sum + wopt(i+1)
      ENDDO
      wopt(1) = 1. - sum 
c
c--------------------------------------------
c     Lagrange multiplier (budget constraint)
c--------------------------------------------
      lambda = dwork(pdlagr + n + neq + nin)    
c
      RETURN
      END
