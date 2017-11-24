c=======================================================================
c
c     subroutine ALLOCSKEW
c
c     Asset allocation with a skewness constraint 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocskew ( n, p, x, w, cov, rho,
     &                       neq, nin, ccst, bcst, cinf, csup, 
     &                       mu, sigma, skew,
     &                       iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of values                               integer
c       p       : portfolio size                                 integer
c       x       : returns values (n*p)                            double
c       w       : benchmark portfolio weights (p)                 double   
c       cov     : covariance matrix (p*p)                         double
c       rho     : mean returns vector (p)                         double
c       neq     : number equality constraints                    integer
c       nin     : number inequality constraints                  integer
c       ccst    : matrix of constraints (p*(neq+nin))             double
c       bcst    : vector initial of constraints (neq+nin)         double
c       cinf    : lower bound (p)                                 double
c       csup    : upper bound (p)                                 double
c       mu      : performance target                              double
c       sigma   : volatility target                               double
c       skew    : skewness delta target (vs. index)               double
c
c     WORKSPACE 
c       iwork   : 16*p + 4*nin + 2*neq + 31                      integer 
c       dwork   : (p+1)*(p+1)*(p+nin+neq+6)+2*nin+neq+2*p+7
c               + (neq + 2*nin + 2*p + 7)*(neq + 2*nin + 2*p + 32)/2
c               + (2*neq + 2*nin + 2*p + 13)*(p + 1)*(p + 2)/2
c               + (p + 1)*(15*p + 31)                             double
c
c     OUTPUT 
c       wopt    : optimal portfolio (p)                           double
c       info    : diagnostic argument                            integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info, neq, nin
      DOUBLE PRECISION mu, sigma, skew
      DOUBLE PRECISION x(*), w(*), cov(*), rho(*), cinf(*), csup(*), 
     &                 ccst(*), bcst(*), wopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwo, pdwo, pisdp, pdsdp, pisdls, pdsdls, pdskew1, 
     &        pdskew3, pdwopt, pdQ, pdCeq, pdbEq, pdCineq, 
     &        pdbLowerIneq, pdbUpperIneq  
      INTEGER i, q, nbEqCst, nbIneqCst
      DOUBLE PRECISION mutest, cstPrecision, minEigenValue
      DOUBLE PRECISION pVolat, bVolat, pSkew, bSkew, pReturn, bReturn, 
     &                 EPS
      PARAMETER ( EPS = 1.E-8 )
c
c     external subroutines
      EXTERNAL opmv, testsdp, utskew1, utskew2, utskew3, sdlsg
c
c-----------------------------------------------------------------------
c
c     initializations 
      info      = 0
      q         = p + 1
      nbEqCst   = neq + 1
      nbIneqCst = nin + p + 1 + 1 + 1
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp : pointer for TESTSDP (p)      
      piwo  = 1
c     piwo  : pointer for OPMV (3*p + 2*nin + neq + 3)
      pisdls = 1
c     pisdls : pointer for SDLSG (2*nbEqCst + 4*nbIneqCst + 12*(p+1) + 5)      
c                               = 2*neq + 4*nin + 4*p + 12*p + 2 + 12 + 12 + 5
c
c     Total size of iwork array = ( 16*p + 4*nin + 2*neq + 31 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdsdp = 1
c     pdsdp : pointer for TESTSDP p*(2*p + 7)           
      pdwo  = 1
c     pdwo  : pointer for OPMV
c            p*( p + neq + nin + 9 ) + 2*neq + 4*nin + 4  
      pdQ          = 1
c     pdQ   : pointer for initial matrix SDLS (p+1)*(p+1)
      pdwopt       = pdQ + (p+1)*(p+1)
c     pdwopt : pointer for optimal matrix SDLS (p+1)*(p+1)     
      pdCeq        = pdwopt + (p+1)*(p+1)
c     pdCeq : pointer for eq. cst. matrix SDLS (neq+1)*(p+1)*(p+1)
      pdbEq        = pdCeq + (neq+1)*(p+1)*(p+1)
c     pdbEq : pointer for eq. cst. vector SDLS (neq+1)      
      pdCineq      = pdbEq + (neq+1)
c     pdCineq : pointer for ineq. cst. matrix SDLS (nin+p+3)*(p+1)*(p+1)       
      pdbLowerIneq = pdCineq + (nin+p+3)*(p+1)*(p+1)
c     pdbLowerIneq : pointer for ineq. cst. vector SDLS (nin+p+3)
      pdbUpperIneq = pdbLowerIneq + (nin+p+3)  
c     pdbUpperIneq : pointer for ineq. cst. vector SDLS (nin+p+3)
      pdskew1      = pdbUpperIneq + (nin+p+3) 
c     pdskew1 : pointer for UTSKEW1 (p*p) 
      pdskew3      = pdbUpperIneq + (nin+p+3) 
c     pdskew3 : pointer for UTSKEW3 p*(15*p + 3) 
      pdsdls       = pdbUpperIneq + (nin+p+3) 
c     pdsdls : pointer for SDLSG           
c                    (mcte+2*mcti)*(mcte+2*mcti+25)/2
c                  + (2*mcte+2*mcti+5)*(p*(p+1)/2)
c                  +  p*(4*p + 27)
c
c     Total size of dwork array = (p + 1)*(p + 1)*(p + nin + neq + 6) + 2*nin + neq + 2*p + 7
c                               + (neq + 2*nin + 2*p + 7)*(neq + 2*nin + 2*p + 32)/2
c                               + (2*neq + 2*nin + 2*p + 13)*(p + 1)*(p + 2)/2
c                               + (p + 1)*(15*p + 31)
c
c-----------------------------------------------------------------------
c
c     test inputs
      IF (sigma .LT. EPS) THEN
        info = -104
        RETURN
      ENDIF
c
c     covariance matrix SDP test
      CALL testsdp (p, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF  
c
c     tests constraints compatibility
      mutest = -1.E+8
      CALL opmv ( p, cov, rho, mutest,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -100
        RETURN
      ENDIF
c
c     test constraint w*rho > mu
      CALL opmv ( p, cov, rho, mu,
     &            neq, nin, ccst, bcst, cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt, info )
      IF (info .EQ. 1001) THEN
        info = -101
        RETURN
      ENDIF
c
c
c     SDLS input matrix
c     Q = | W   w |
c         | w'  1 |
      CALL utskew1(p, w, iwork(pisdp), dwork(pdskew1), dwork(pdQ), info)
      IF (info .LT. 0) RETURN
c  
c     SDLS contraints: equalities      
      CALL utskew2(p, neq, ccst, bcst,
     &             nbEqCst, dwork(pdCeq), dwork(pdbEq), info)
      IF (info .LT. 0) RETURN
c
c     SDLS contraints: inequalities 
      CALL utskew3(n, p, x, w, rho, cov,
     &             neq, nin, ccst, bcst, cinf, csup,
     &             mu, sigma, skew, 
     &             iwork(pisdp), dwork(pdskew3), 
     &             nbIneqCst, dwork(pdCineq), 
     &             dwork(pdbLowerIneq), dwork(pdbUpperIneq), info)
      IF (info .LT. 0) RETURN
c
c     optimization with SDLS
      cstPrecision  = 1.E-15
      minEigenValue = 1.E-12
      CALL sdlsg ( q, dwork(pdQ), nbEqCst, dwork(pdCeq), dwork(pdbEq),
     &             nbIneqCst, dwork(pdCineq), 
     &             dwork(pdbLowerIneq), dwork(pdbUpperIneq),
     &             cstPrecision, minEigenValue,
     &             iwork(pisdls), dwork(pdsdls), dwork(pdwopt), info )
      IF (info .NE. 0) RETURN
c
c     optimal portfolio
      DO i = 1,p
        wopt(i) = dwork(pdwopt + q*q - q - 1 + i)
      ENDDO
c
c     ex-ante skewness (benchmark)
      CALL EXARSK ( n, p, w, x, cov, dwork(pdsdls), bSkew, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante skewness (portfolio)
      CALL EXARSK ( n, p, wopt, x, cov, dwork(pdsdls), pSkew, info)
      IF (info .LT. 0) RETURN
c      
      CALL EXARET ( p, rho, w, bReturn, info ) 
      IF (info .LT. 0) RETURN
c      
c     ex-ante mean return
      CALL EXARET ( p, rho, wopt, pReturn, info ) 
      IF (info .LT. 0) RETURN
c
c     ex-ante volatility (optimal portfolio)      
      CALL EXARVO (p, cov, w, bVolat, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante volatility (optimal portfolio)      
      CALL EXARVO (p, cov, wopt, pVolat, info)
      IF (info .LT. 0) RETURN
c      
      IF (pSkew .LT. (bSkew+skew)) THEN 
        info = 120
        RETURN
      ENDIF
      IF (pReturn .LT. (mu - EPS)) THEN
        info = 121
        RETURN
      ENDIF    
      IF (pVolat .GT. (sigma + EPS)) THEN 
        info = 122
        RETURN
      ENDIF
c
      RETURN
      END
