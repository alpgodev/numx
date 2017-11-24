c========================================================================
c
c     subroutine OPVAR
c
c     Computing optimization for ALLOCVAR - Value-at-Risk constraint
c
c-----------------------------------------------------------------------
      SUBROUTINE opvar ( n, cov, rho, neq, nin, ccst, bcst,
     &                   cinf, csup, vrimax, conflp, maxdic, epsdic,
     &                   iwork, dwork, wopt, vriopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of assets                           integer
c            cov    : covariance matrix (n*n)                     double
c            rho    : mean returns vector (n)                     double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix of constraints (n*(neq+nin))         double
c            bcst   : vector of constraints  vector(neq+nin)      double
c            cinf   : lower bounds (n)                            double
c            csup   : upper bounds (n)                            double
c            vrimax : Value-at-Risk constraint                    double
c            conflp : confidence level parameter ( 0<conflp<1 )   double
c            maxdic : maximum of dichotomy iterations            integer
c            epsdic : precision of dichotomy                      double
c
c     WORKSPACE 
c            iwork  : 3*n + neq + 2*nin + 3                      integer 
c            dwork  : n*(n + neq + nin + 10) + 2*neq + 4*nin + 4  double
c
c     OUTPUT 
c            wopt   : optimal portfolio (n)                       double
c            vriopt : optimal Value-at-Risk                       double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            EVBOR  : computing the minimum and maximum of a vector
c            OPMV   : computing optimization for ALLOCMV
c            WZERO  : put zero at the values of a vector less than eps
c                     and increase the greatest element with residual
c            YV     : copy a vector in a vector
c            XV     : computing scalar product of two vectors = V1'*V2
c            OVTMCV : computing V'*M*V = scalar
c                    ( M vectorized square matrix(n*n), V vector(n) )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, maxdic, info
      DOUBLE PRECISION vrimax, conflp, epsdic, vriopt
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER niter, infmax, infqua, piwo, pdwopt, pdwo
      DOUBLE PRECISION p, q, zc, rhopt, optmean, optvar, volat, vrisav,
     &                 rhomin, rhomax, epsvri, xvri, eps, ornumber
      PARAMETER ( eps = 1.0E-15 )
c
c     external subroutines
      EXTERNAL EVBOR, opmv, WZERO, YV, XV, OVTMCV, dinvnr
      DOUBLE PRECISION dinvnr
c      
c     intrinsic functions
      INTRINSIC sqrt, min, max, abs         
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwo = 1
c     piwo  : pointer for OPMV ( 3*n + neq + 2*nin + 3 )
c
c     Total size of iwork array = ( 3*n + neq + 2*nin + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwopt = 1
c     pdwopt : pointer for wopt save vector(n)
      pdwo   = pdwopt + ( n )
c     pdwo   : pointer for internal workspace of OPMV
c              ( n*( n + neq + nin + 9 ) + 2*neq + 4*nin + 4  )
c
c     Total so size of dwork array = n
c                                  + n*(n+neq+nin+9) 
c                                  + 2*neq + 4*nin + 4
c         = n*(n + neq + nin + 10) + 2*neq + 4*nin + 4
c
c-----------------------------------------------------------------------
c
c     rdt min/max (rhomin and rhomax)
      CALL EVBOR ( n, rho, rhomin, rhomax )
c
c     research optimal performance subject to Value-at-Risk constraint 
c     (upper bound VaR) by dichotomy
c
c     initializations
      vrisav = 0.
      rhopt  = rhomax
      epsvri = 2.*epsdic
      vriopt = -1.0
      infmax = -107
      infqua = 0
      niter  = 0
      xvri   = -1.0
c
c     or number
      ornumber = 2.0/(1.0 + SQRT(5.0))
c
c     loop dichotomy search
      DO WHILE ( ( (niter .LE. maxdic).AND.(epsvri .GT. epsdic) ) .AND.
     &           ( (vrimax .GT. xvri) .OR. (niter .NE. 1) ) )
         niter = niter + 1
c
c        Mean-Variance allocation (Markowitz)
         CALL opmv ( n, cov, rho, rhopt,
     &               neq, nin, ccst, bcst, cinf, csup,
     &               iwork(piwo), dwork(pdwo), wopt, info )
c
c        if wopt(i)< eps then wopt(i) = 0
c        the residual increase max[wopt(i)]
         CALL WZERO ( n, wopt, eps )
c
c        error management 
         IF (info .NE. 0) THEN
            IF (info .NE. 1001) THEN
                RETURN
            ELSE
                infqua = 1001
            ENDIF
         ELSE
            infqua = 0
            CALL YV ( n, wopt, dwork(pdwopt) )
         ENDIF
c
c        optimal portfolio ex-ante variance = w'*cov*w
         CALL OVTMCV ( n, cov, wopt, optvar )
c
c        optimal portfolio ex-ante volatility = sqrt(w'*cov*w)
         volat = sqrt(optvar)
c
c        optimal ex-ante mean return = w'*rmean
         CALL XV ( n, wopt, rho, optmean )
c
c        inverse of the N(0,1) cumulative distribution
c        call dinvnr(p,q)
         p  = conflp
         q  = 1. - p
         zc = dinvnr(p,q)
c
c     write p, q, zc
c      open(unit=1,file='zc.txt',status='unknown')
c	write(1,*) p
c	write(1,*) q
c	write(1,*) zc
c	close(unit=1)
c
c        ex-ante Value-at-Risk (upper bound)
c         xvri = - volat * sqrt( conflp/(1.-conflp) )
c         xvri = - volat * sqrt( conflp/(1.-conflp) ) + emret
c     
c        ex-ante Normal Value-at-Risk 
         xvri = - volat * zc
c         xvri = - volat * zc + emret       
c
c     write VaR
c      open(unit=1,file='VaR.txt',status='unknown')
c	write(1,*)  xvri
c	close(unit=1)
c
c        impose that VaR in [-100%, 100%]
         xvri = MIN(xvri, 1.)
         xvri = MAX(xvri, -1.)
c
c        save the last best result
         IF (infqua .EQ. 0) THEN
            vrisav = xvri
         ELSE
            IF (niter .EQ. 1) THEN
               vrisav = xvri
               CALL YV ( n, wopt, dwork(pdwopt) )
               infqua = 0
               info   = 0
            ENDIF
         ENDIF
c
c        stop test
         epsvri = abs(xvri - vriopt)
         vriopt = xvri
c
c        dichotomy (objective return) 
         IF (vrimax.lt.vriopt) THEN
            rhomin = rhopt
            infmax = 0
         ELSE
            rhomax = rhopt
         ENDIF
         rhopt = ( rhomax + rhomin ) / 2.
c          rhopt = ornumber * ( rhomax + rhomin )
      ENDDO
c
c     output management
      IF (niter.gt.maxdic) info = 106
      IF (infmax.ne.0) info = -infmax
      IF (infqua.ne.0) THEN
         info = infqua
         CALL YV ( n, dwork(pdwopt), wopt )
         vriopt = vrisav
      ENDIF
c
      RETURN
      END
