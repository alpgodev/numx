c=======================================================================
c
c     subroutine OPIT                                        
c     
c     Optimization utility for ALLOCIT - Index Tracking Allocation
c
c-----------------------------------------------------------------------
      SUBROUTINE opit ( n, cov, covb, rho, rhob, delta,
     &                  neq, nin, ccst, bcst, cinf, csup,
     &                  iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : size of portfolio                          integer
c            cov    : covariance matrix (n*n)                     double
c            covb   : covariance asset/benchmark (n)              double
c            rho    : asset mean returns vector (n)               double
c            rhob   : bechmark mean returns                       double
c            delta  : relative target return                      double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix of constraints (n*(neq + nin))       double
c            bcst   : vector of constraints (neq + nin)           double
c            cinf   : assets inferior limit (n)                   double
c            csup   : assets superior limit (n)                   double
c
c     WORKSPACE 
c
c
c            iwork  : 3*n+2*nin+neq+6                            integer
c            dwork  : (n+1)*(neq+nin+2*(n+1)+11)+neq+3*nin        double
c
c     OUTPUT 
c            wopt   : optimal portfolio (n)                       double
c            info   : diagnostic argument                        integer
c
c     CALL   
c       CTIT     : build constraints
c       INITFEAS :
c       QP       : quadratic solver
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION delta, rhob
      DOUBLE PRECISION cov(*), covb(*), rho(*), cinf(*), csup(*), 
     &                 wopt(*), ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER neqtot, nintot, piw, pdlin, pdw, pdlagr, pdbcs, pdccs
c
c     external subroutines
      EXTERNAL ctit, initfeas, qp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw  = 1
c     piw    : pointer for initfeas 
c               needs 3*n1 + 2*nin + neq + 1
c               so 3*n+2*nin+neq+6
c      
c     piw  : pointer for QP internal workspace, who needs
c            ( 3*n + 2*nintot + neqtot + 1 )
c             with neqtot = neq
c                  nintot = nin + 1
c
c     Total size of iwork array = ( 3*n + 2*nin + neq + 6 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n ) more
      pdlagr    = pdlin + ( n )
c     pdlagr : pointer for Lagrange multipliers vector
c              ( n + nintot + neqtot ),
c           so ( n + nin + 1 + neq ) more
      pdbcs  = pdlagr +  ( n + nin + 1 + neq )
c     pdbcs  : pointer for the constraints vector ( nintot + neqtot ),
c              so ( nin + 1 + neq ) more
      pdccs  = pdbcs +  ( nin + 1 + neq )
c     pdccs  : pointer for the constraints matrix
c              n*( nintot + neqtot ),
c           so n*( nin + 1 + neq ) more
      pdw    = pdccs + ( n*( nin + 1 + neq ) )
c     pdw    : pointer for initfeas 
c               needs (n+1)*(neq+nin+2*(n+1)+11)+neq+3*nin  
c     pdw    : pointer for QP internal workspace who needs
c              (n*n + 6*n + 2*nintot)
c               so pdw needs (n+1)*(neq+nin+2*(n+1)+11)+neq+3*nin  
c
c     Total size of dwork array = n
c                               + n + nin + 1 + neq 
c                               + nin + 1 + neq 
c                               + n*( nin + 1 + neq ) 
c                               + n*(neq+nin+2*n+15)+2*neq+4*nin+13 
c                               = n*(2*n + 2*nin + 2*neq + 18) + 4*neq + 6*nin + 15
c
c      = n*(2*n + 2*nin + 2*neq + 18) + 4*neq + 6*nin + 15
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector and linear part
      CALL ctit ( n, rhob, rho, covb, delta, neq, nin, ccst, bcst,
     &            neqtot, nintot, dwork(pdccs), dwork(pdbcs),
     &            dwork(pdlin) )
c
c     problem initialization and feasibility check
      CALL initfeas ( n, cinf, csup, neqtot, nintot, dwork(pdccs),
     &                dwork(pdbcs), iwork(piw), dwork(pdw),
     &                wopt, info)
      IF (info .NE. 0) THEN
        info = -100
        RETURN
      ENDIF
c
c     quadratic solver   
      CALL qp ( n, cov, dwork(pdlin), neq, nintot,
     &          dwork(pdccs), dwork(pdbcs), cinf, csup,
     &          iwork(piw), dwork(pdw), dwork(pdlagr), wopt,
     &          info )
      RETURN
      END
