c=======================================================================
c
c     subroutine OPVOLRFR                                      
c     
c     Optimization utility for ALLOCVOLRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE opvolrfr ( n, cov, rho, rfr, mu, neq, nin, ccst, bcst,
     &                      cinf, csup, cinfrfr, csuprfr, 
     &                      iwork, dwork, wopt, lagr, info )
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
c       dwork  : n*(n+neq+nin+10) + neq  + 3*nin + 9               double

c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       lagr   : Lagrange (n + neq + 3 + nin)                     double
c       info   : = 0 successful exit                             integer
c
c     CALL   
c       CTMVRFR : utility function (cf. modules/UTALLOC.F)
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
     &                 ccst(*), bcst(*), lagr(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER neqtot, nintot, piw, pdlin, pdw, pdbcs, pdccs
c
c     external subroutines
      EXTERNAL ctmvrfr, qp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw  = 1
c     piw  : pointer for QP internal workspace, who needs
c            ( 3*n + 2*nintot + neqtot + 1 )
c             with neqtot = neq
c                  nintot = nin + 3
c
c     Total size of iwork array = 3*n + 2*nin + neq + 7 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n ) more
      pdbcs  = pdlin +  ( n )
c     pdbcs  : pointer for the constraints vector ( nintot + neqtot ),
c              so ( nin + 3 + neq ) more
      pdccs  = pdbcs +  ( nin + 3 + neq )
c     pdccs  : pointer for the constraints matrix
c              n*( nintot + neqtot ),
c           so n*( nin + 3 + neq ) more
      pdw    = pdccs + ( n*( nin + 3 + neq ) )
c     pdw    : pointer for QP internal workspace who needs
c              (n*n + 6*n + 2*nintot)
c
c     Total size of dwork array = n
c                               + nin + 3 + neq 
c                               + n*( nin + 3 + neq ) 
c                               + n*n + 6*n + 2*nin + 6    
c                               = n*(n + neq + nin + 10) + neq  + 3*nin + 9  
c
c      = n*(n + neq + nin + 10) + neq  + 3*nin + 9
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector, and linear part
      CALL ctmvrfr ( n, rho, rfr, mu, neq, nin, ccst, bcst,
     &               cinfrfr, csuprfr, neqtot, nintot,
     &               dwork(pdccs), dwork(pdbcs), dwork(pdlin))
c
c     quadratic solver
      CALL qp ( n, cov, dwork(pdlin),
     &          neqtot, nintot, dwork(pdccs), dwork(pdbcs),
     &          cinf, csup, iwork(piw), dwork(pdw),
     &          lagr, wopt, info )
c
      RETURN
      END
