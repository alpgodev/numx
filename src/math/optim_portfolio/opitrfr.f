c=======================================================================
c
c     subroutine OPITRFR
c     
c     Optimization utility for ALLOCITRFR - Index Tracking Allocation
c
c-----------------------------------------------------------------------
      SUBROUTINE opitrfr ( n, cov, covb, rho, rhob, rfr, delta,
     &                     neq, nin, ccst, bcst, 
     &                     cinf, csup, cinfrfr, csuprfr,
     &                     iwork, dwork, wopt, lagr, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : size of portfolio                          integer
c            cov    : covariance matrix (n*n)                     double
c            covb   : covariance asset/benchmark (n)              double
c            rho    : asset mean returns vector (n)               double
c            rhob   : bechmark mean returns                       double
c            rfr    : risk-free rate                              double
c            delta  : relative target return                      double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix of constraints (n*(neq + nin))       double
c            bcst   : vector of constraints (neq + nin)           double
c            cinf   : assets inferior limit (n)                   double
c            csup   : assets superior limit (n)                   double
c       cinfrfr: risk-free rate lower bound                       double
c       csuprfr: risk-free rate upper bound                       double
c
c     WORKSPACE 
c            iwork  : 3*n + 2*nin + neq + 7                      integer 
c            dwork  : n*(2*n+2*nin+2*neq+22)+3*neq+5*nin+28       double

c     OUTPUT 
c            wopt   : optimal portfolio (n)                       double
c            lagr   : Lagrange (n + neq + 3 + nin)                double
c            info   : diagnostic argument                        integer
c
c     CALL   
c       CTITRFR  : build constraints
c       INITFEAS :
c       QP       : quadratic solver
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION delta, rhob, rfr, cinfrfr, csuprfr
      DOUBLE PRECISION cov(*), covb(*), rho(*), cinf(*), csup(*), 
     &                 wopt(*), lagr(*), ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, ntot, neqtot, nintot
      INTEGER piw, pdlin, pdw, pdb, pdc
      DOUBLE PRECISION ZERO, EPS
      PARAMETER (ZERO = 0.D0, EPS = 1.E-8)
c
c     external subroutines
      EXTERNAL ctitrfr, initfeas, qp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw  = 1
c     piw  : pointer for QUAPRO   3*n + 2*nin + neq + 7
c            pointer for INITFEAS 3*n + 2*nin + neq + 7
c
c     Total size of iwork array = 3*n + 2*nin + neq + 7
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n ) more
      pdb    = pdlin +  ( n )
c     pdbcs  : pointer for the constraints vector ( neq + nin + 3 )
      pdc    = pdb +  ( neq + nin + 3 )
c     pdccs  : pointer for the constraints matrix n*( neq + nin + 3 )
      pdw    = pdc + ( n*( neq + nin + 3 ) )
c     pdw    : pointer for QUAPRO   (n*n + 6*n + 2*nin + 6)
c                      for INITFEAS  n*(2*n+neq+nin+18)+4*nin+2*neq+25
c
c     Total size of dwork array = n
c                               + nin + neq + 3 
c                               + n*( nin + neq + 3 ) 
c                               + n*(neq+nin+2*n+18)+2*neq+4*nin+25
c                               = n*(2*n+2*nin+2*neq+22)+3*neq+5*nin+28
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector and linear part
      ntot   = neq + nin  ! number total of input constraint
      neqtot = neq        ! number of equalities constraints 
      nintot = nin + 3    ! number of inequalities constraints
c
c     construction of the linear part
      DO i = 1,n
         dwork(pdlin + i - 1) = -covb(i)
      ENDDO
c
c     case ncst <> 0 (neq=0, nin=0)
      IF (ntot .NE. 0) THEN
        CALL YV ( n*ntot, ccst, dwork(pdc) )  ! constraint matrix
        CALL YV ( ntot, bcst, dwork(pdb) )    ! constraint vector
      ENDIF
c 
      DO i = 1,n
        dwork(pdc + ntot*n + i - 1)     = rfr - rho(i)
        dwork(pdc + (ntot+1)*n + i - 1) =  1.0
        dwork(pdc + (ntot+2)*n + i - 1) = -1.0
      ENDDO
      dwork(pdb + ntot)     = rfr - delta - rhob
      dwork(pdb + ntot + 1) = 1.0 - cinfrfr + EPS
      dwork(pdb + ntot + 2) = csuprfr - 1.0 + EPS
c
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neq, nintot,
     &                dwork(pdc), dwork(pdb),
     &                iwork(piw), dwork(pdw),
     &                wopt, info)
      IF (info .EQ. 1111) THEN
        info = -100
        RETURN
      ENDIF
      IF (info .NE. 0) RETURN
c     
c     quadratic solver
      CALL qp ( n, cov, dwork(pdlin), neq, nintot,
     &          dwork(pdc), dwork(pdb), cinf, csup,
     &          iwork(piw), dwork(pdw), lagr, wopt,
     &          info )
      RETURN
      END
