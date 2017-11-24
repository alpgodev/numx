c=======================================================================
c
c     subroutine ALLOCVOL                            
c
c     This function implements Volatility-constrained optimization
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
      SUBROUTINE allocvol ( n, cov, rho,
     &                      neq, nin, ccst, bcst, cinf, csup, 
     &                      volb, volp,
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
c       volb   : volatility budget (>0)                           double 
c       volp   : precision of the volatility cst. (>0)            double
c
c     WORKSPACE 
c       iwork  : 3*n + 2*nin + neq + 4                          integer 
c       dwork  : n*(2*n+neq+nin+9) + 4*nin + 2*neq + 5          double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
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
      DOUBLE PRECISION epskat, volb, volp, lambda
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(n,*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, niter, maxiter,
     &        pisdp, piwo, pdsdp, pdwo, pdlagr
      DOUBLE PRECISION mu, EPS, rhomax, rhomin, rhomean, tmin, tmax, t,
     &                 error, mincinf, var, volmin, volopt, epsilon
      PARAMETER (EPS = 1.E-50)
c
c     external subroutines
      EXTERNAL EVMAX, EVMIN, YM, opvol, testsdp, OVTMCV
c     
c     intrinsic functions
      INTRINSIC MAX, ABS, SQRT
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
c             needs( 3*n + 2*nin + neq + 3 )
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP (n*(2*n + 7))              
      pdwo  = 1
c     pdwo  : pointer for internal workspaces of OPVOL
c             needs( n*( n + neq + nin + 8 ) 
c                    + neq + 3*nin + 3  )
      pdlagr = pdwo + ( n*(n+neq+nin+8) + neq + 3*nin + 3 )
c     pdlagr : pointer for Lagr of OPVOL
c              needs ( n + neq + nin + 1 )      
c
c     Total size of dwork array = n*( n + neq + nin + 8 )
c                               + neq + 3*nin + 3 
c                               + n + neq + nin + 1
c     = n*( 2*n + neq + nin + 9 ) + 4*nin + 2*neq + 4
c
c
c-----------------------------------------------------------------------
c      open(unit=1,file='ALLOCVOL.txt',status='unknown')
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
c     test if budget volatility >= min. volatility
      t = -1.E+15
c      CALL OPMV ( n, cov, rho, t,
c     &            neq, nin, ccst, bcst, cinf, csup,
c     &            iwork(piwo), dwork(pdwo), wopt, info )
c       optimization (cf. OPVOL.F)     
       CALL opvol ( n, cov, rho, t,
     &              neq, nin, ccst, bcst, cinf, csup,
     &              iwork(piwo), dwork(pdwo), wopt, dwork(pdlagr), info)
c
c     if (w'*Q*w) < volb2 -> volatility budget const. too small !
      CALL OVTMCV ( n, cov, wopt, var)
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
        IF (ABS(rho(i)) .LT. EPS) THEN
            info = -110
            RETURN
        ENDIF
      ENDDO
      
c------------------------------------
c     initialization of tmin and tmax
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
        CALL EVMIN(n, rho, rhomin)
        tmax = rhomax
        tmin = rhomin
c
c     case long/short positions        
      ELSE
        rhomax = ABS(rho(1))
        DO i = 2,n
            rhomax = MAX(rhomax, ABS(rho(i))) 
        ENDDO
        tmax = 0.0
        DO i = 1,n
            tmax = tmax + rhomax * (MAX(ABS(cinf(i)), ABS(csup(i))))  
        ENDDO
        tmin = -tmax
      ENDIF
c
c---------------------------------
c     dichotomy
c---------------------------------
c
c     myzero := max(1.E-12, mean(rho))
      CALL MV ( n, rho, rhomean )
      epsilon = MAX(1.E-12, rhomean*1.E-8)
c
c     dichotomy loop
      error = 1.E+15
      niter = 0
c      
c     test if the volatility budget constraint is too high 
c     (budget constraint not attainable)
      IF ( ABS(tmax - tmin) .LT. epsilon ) THEN
        info = 112
        RETURN
      ENDIF 
c      
      DO WHILE ((niter .LT. maxiter) .AND. (error .GT. volp))
c 
        niter = niter + 1
        t = (tmin + tmax)/2
c
c       optimization (cf. opvol.f)
        CALL opvol ( n, cov, rho, t,
     &          neq, nin, ccst, bcst, cinf, csup,
     &          iwork(piwo), dwork(pdwo), wopt, dwork(pdlagr), info )
c
c       sqrt[w'*Q*w] (portfolio volatility)
        CALL OVTMCV ( n, cov, wopt, var)      
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
c
c     maximum iterations      
      IF (niter .EQ. maxiter) THEN
        info = 111
      ENDIF
c      
c      close(unit=1)
c
c--------------------------------------------
c     Lagrange multiplier (budget constraint)
c--------------------------------------------
c
      mu = dwork(pdlagr + n + neq + nin)
      IF (ABS(mu) .LT. EPS) THEN
        lambda = 1./EPS  
      ELSE  
        lambda = 1./mu
      ENDIF
c
      RETURN
      END
