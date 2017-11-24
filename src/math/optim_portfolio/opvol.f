c=======================================================================
c
c     subroutine OPVOL                                      
c     
c     Optimization utility for ALLOCVOL
c
c-----------------------------------------------------------------------
      SUBROUTINE opvol ( n, cov, rho, mu, neq, nin, ccst, bcst,
     &                  cinf, csup,
     &                  iwork, dwork, wopt, lagr, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : size of portfolio                          integer
c            cov    : covariance matrix (n*n)                     double
c            rho    : mean returns vector (n)                     double
c            mu     : performance target                          double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix of constraints (n*(neq + nin))       double
c            bcst   : vector of constraints (neq + nin)           double
c            cinf   : assets inferior limit (n)                   double
c            csup   : assets superior limit (n)                   double
c
c     WORKSPACE 
c            iwork  : 3*n + 2*nin + neq + 3                      integer 
c            dwork  : n*(n+neq+nin+8) + neq  + 3*nin + 3          double

c     OUTPUT 
c            wopt   : optimal portfolio (n)                       double
c            lagr   : Lagrange (n + neq + 1 + nin)                double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c            CTMV : utility function (cf. utalloc.f)
c            QP   : quadratic solver (vectorized version)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION mu
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), wopt(*), 
     &                 ccst(n,*), bcst(*), lagr(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER neqtot, nintot, piw, pdlin, pdw, pdbcs, pdccs
c
c     external subroutines
      EXTERNAL ctmv, qp
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
c                  nintot = nin + 1
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n ) more
      pdbcs  = pdlin +  ( n )
c     pdbcs  : pointer for the constraints vector ( nintot + neqtot ),
c              so ( nin + 1 + neq ) more
      pdccs  = pdbcs +  ( nin + 1 + neq )
c     pdccs  : pointer for the constraints matrix
c              n*( nintot + neqtot ),
c           so n*( nin + 1 + neq ) more
      pdw    = pdccs + ( n*( nin + 1 + neq ) )
c     pdw    : pointer for QP internal workspace who needs
c              (n*n + 6*n + 2*nintot)
c
c     Total size of dwork array = n
c                               + nin + 1 + neq 
c                               + n*( nin + 1 + neq ) 
c                               + n*(n + 6) + 2*(nin + 1)    
c                               = n*(n + neq + nin + 8) + neq  + 3*nin + 3  
c
c      = n*(n + neq + nin + 8) + neq  + 3*nin + 3
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector, and linear part
      CALL ctmv ( n, rho, mu, neq, nin, ccst, bcst, neqtot, nintot,
     &            dwork(pdccs), dwork(pdbcs), dwork(pdlin))
c
c     quadratic solver
      CALL qp ( n, cov, dwork(pdlin), dwork(pdccs),
     &          dwork(pdbcs), cinf, csup, neqtot, nintot,
     &          iwork(piw), dwork(pdw), lagr, wopt,
     &          info )
c
      RETURN
      END
