c=======================================================================
c
c     subroutine OPMVRFR      
c     
c     Optimization utility for ALLOCMVRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE opmvrfr ( n, cov, rho, rfr, mu, neq, nin, ccst, bcst,
     &                     cinf, csup, cinfrfr, csuprfr, 
     &                     iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size of portfolio                               integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : expected mean returns vector (n)                 double
c       rfr    : risk-free rate                                   double
c       mu     : performance target                               double
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       ccst   : matrix of constraints (n*(neq + nin))            double
c       bcst   : vector of constraints (neq + nin)                double
c       cinf   : assets inferior limit (n)                        double
c       csup   : assets superior limit (n)                        double
c       cinfrfr: risk-free rate lower bound                       double
c       csuprfr: risk-free rate upper bound                       double
c
c     WORKSPACE 
c       iwork  : 3*n + 2*nin + neq + 7                           integer 
c       dwork  : n*(n+neq+nin+11)+2*neq+4*nin+12                  double

c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       info   : = 0 successful exit                             integer
c
c     CALL   
c       QP      : quadratic solver (vectorized version)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION mu, rfr, cinfrfr, csuprfr
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), wopt(*), 
     &                 ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, ntot, neqtot, nintot
      INTEGER piw, pdlin, pdw, pdb, pdc, pdlagr
      DOUBLE PRECISION ZERO, EPS
      PARAMETER (ZERO = 0.D0, EPS = 1.E-8)
c
c     external subroutines
      EXTERNAL qp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw  = 1
c     piw  : pointer for QUAPRO ( 3*n + neq + 2*nin + 7 )
c
c     Total size of iwork array = 3*n + 2*nin + neq + 7 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n ) more
      pdlagr = pdlin + ( n )
c     pdlagr : pointer for Lagrange multipliers ( n + neq + nin + 3 )
      pdb    = pdlagr + ( n + neq + nin + 3 )
c     pdbcs  : pointer for the constraints vector ( neq + nin + 3 )
      pdc    = pdb +  ( neq + nin + 3 )
c     pdccs  : pointer for the constraints matrix n*( neq + nin + 3 )
      pdw    = pdc + ( n*( neq + nin + 3 ) )
c     pdw    : pointer for QUAPRO  (n*n + 6*n + 2*nin + 6)
c
c     Total size of dwork array = n
c                               + n + neq + nin + 3
c                               +     neq + nin + 3 
c                               + n*( neq + nin + 3 ) 
c                               + n*n + 6*n + 2*nin + 6    
c
c     Total size of dwork array = n*(n + neq + nin + 11) + 2*neq  + 4*nin + 1é  
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector and linear part
      ntot   = neq + nin  ! number total of input constraint
      neqtot = neq
      nintot = nin + 3
c
c     construction of the linear part
      CALL IVX (n, dwork(pdlin), ZERO)
c
c     case ncst <> 0 (neq=0, nin=0)
      IF (ntot .NE. 0) THEN
        CALL YV ( n*ntot, ccst, dwork(pdc) )  ! constraint matrix
        CALL YV ( ntot, bcst, dwork(pdb) )    ! constraint vector
      ENDIF
      DO i = 1,n
        dwork(pdc + ntot*n + i - 1)     =  rfr - rho(i)
        dwork(pdc + (ntot+1)*n + i - 1) =  1.0
        dwork(pdc + (ntot+2)*n + i - 1) = -1.0
      ENDDO
      dwork(pdb + ntot)     = rfr - mu
      dwork(pdb + ntot + 1) = 1.0 - cinfrfr + EPS
      dwork(pdb + ntot + 2) = csuprfr - 1.0 + EPS
c
c     quadratic solver
      CALL qp ( n, cov, dwork(pdlin), neqtot, nintot,
     &          dwork(pdc), dwork(pdb), cinf, csup,
     &          iwork(piw), dwork(pdw), dwork(pdlagr), wopt,
     &          info )
      RETURN
      END
