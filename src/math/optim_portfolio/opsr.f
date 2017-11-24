c=======================================================================
c
c     subroutine OPSR
c
c     Optimization utility for ALLOCSR - Sharpe Ratio maximization
c
c-----------------------------------------------------------------------
      SUBROUTINE opsr ( n, cov, rho, rfr,
     &                  neq, nin, ccst, bcst, cinf, csup,
     &                  iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size of portfolio                               integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : mean returns vector (n)                          double
c       rfr    : risk-free rate                                   double
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       ccst   : matrix of constraints (n*(neq + nin))            double
c       bcst   : vector of constraints (neq + nin)                double
c       cinf   : assets inferior limit (n)                        double
c       csup   : assets superior limit (n)                        double
c
c     WORKSPACE 
c       iwork  : 5*n + 2*nin + neq + 4                           integer 
c       dwork  : n*(4*n + nin + neq + 25) + 5*nin + 3*neq + 15    double
c
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       info   : = 0 successful exit                             integer
c
c     CALL   
c       CTSR : utility function
c       QP   : quadratic solver (vectorized version)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION rfr
      DOUBLE PRECISION cov(*), rmean(*), cinf(*), csup(*), wopt(*), 
     &                 ccst(n,*), bcst(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, neqtot, nintot, piw, pdlin, pdw, pdlagr, pdbcs, pdccs, 
     &        pdcov, pdinf, pdsup, j, k
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 1.E-6)
c
c     external subroutines
      EXTERNAL ctsr, qp
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
c            ( 3*(n+1) + 2*nintot + neqtot + 1 )
c             with neqtot = neq + 1
c                  nintot = nin + 2*n
c
c     Total size of iwork array = ( 5*n + 2*nin + neq + 4 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------  
      pdlin  = 1
c     pdlin  : pointer for linear part vector, so ( n + 1 ) more
      pdlagr = pdlin + ( n + 1 )
c     pdlagr : pointer for Lagrange multipliers vector
c              ( n + 1 + nintot + neqtot ),
c           so ( n + 1 + nin + 1 + neq + 2*n) more
      pdbcs  = pdlagr +  ( 3*n + nin + neq + 2)
c     pdbcs  : pointer for the constraints vector ( nintot + neqtot ),
c              so ( nin + 1 + neq + 2*n) more
      pdccs  = pdbcs +  ( nin + neq + 1 + 2*n)
c     pdccs  : pointer for the constraints matrix
c              ( n + 1 )*( nintot + neqtot ),
c           so ( n + 1) *( nin + 1 + neq + 2*n) more
      pdw    = pdccs + ( (n + 1)*( nin + neq + 1 + 2*n) )
c     pdw    : pointer for QP internal workspace who needs
c              ((n+1)*(n+1) + 6*(n+1) + 2*(nin + 2*n))
      pdcov  = pdw + ((n+1)*(n+1) + 6*(n+1) + 2*(nin + 2*n))
c     pdcov  : pointer for temporary cov. matrix (n+1)*(n+1)
      pdinf  = pdcov + ((n+1)*(n+1))
c     pdinf  : pointer for temporary lower bounds (n+1)
      pdsup  = pdinf + (n + 1)
c     pdsup  : pointer for temporary upper bounds (n+1)      
c
c     Total size of dwork array = n + 1
c                               + 3*n + nin + neq + 2
c                               + 2*n + nin + neq + 1 
c                               + (n + 1)*( nin + neq + 1 + 2*n )
c                               + (n + 1)*(n + 1) + 6*(n+1) + 2*(nin + 2*n)
c                               + (n + 1)*(n + 1)
c                               + (n + 1) 
c                               + (n + 1)
c                               = n*(4*n + nin + neq + 25) + 5*nin + 3*neq + 15
c
c-----------------------------------------------------------------------
c
c     construction of the constraints matrix and vector, and linear part
      CALL ctsr ( n, rmean, rfr, neq, nin, ccst, bcst, cinf, csup,
     &            neqtot, nintot, dwork(pdccs), dwork(pdbcs), 
     &            dwork(pdlin))
c
c     new cov. matrix cov = | cov  0  |  
c                           |  0  eps | 
      DO i = 1,((n+1)*(n+1))
         dwork(pdcov + i) = 0.
      ENDDO   
      DO i = 1,n
        DO j = 1,n
            k = i + n*(j-1)
            dwork(pdcov + i + (n+1)*(j-1)) = cov(k)
        ENDDO
      ENDDO
      dwork(pdcov + (n+1)*(n+1)) = eps
c
c     lower/upper bounds of the new problem
      DO i = 1,n
        dwork(pdinf + i) =  0.
        dwork(pdsup + i) =  100000.
      ENDDO
      dwork(pdinf + n + 1) = 0.
      dwork(pdsup + n + 1) = 100000.
c
c     write input data
c     write C and b
c      ioecr = 10
c      i = n + 1
c      open(unit=1,file='n.txt',status='unknown')
c	write(1,*) neqtot
c	write(1,*) nintot
c	write(1,*) n
c	close(unit=1)
c
c      open(ioecr,status='unknown',file='cov.txt')
c      call IMPRM ( ioecr, 'cov :', i, i, dwork(pdcov+1) )
c      open(ioecr,status='unknown',file='lin.txt')
c      call IMPRV ( ioecr, 'linear part :', i, dwork(pdlin) )
c      
c      open(ioecr,status='unknown',file='C.txt')
c      call IMPRM ( ioecr, 'C :', n+1, neqtot + nintot, dwork(pdccs) )
c      open(ioecr,status='unknown',file='b.txt')
c      call IMPRV ( ioecr, 'b :', neqtot + nintot, dwork(pdbcs) )
c      
c      open(ioecr,status='unknown',file='cinf.txt')
c      call IMPRV ( ioecr, 'Cinf :', i, dwork(pdinf+1) )
c      open(ioecr,status='unknown',file='csup.txt')
c      call IMPRV ( ioecr, 'Csup :', i, dwork(pdsup+1) )
c
c     quadratic solver
      CALL qp ( (n+1), dwork(pdcov+1), dwork(pdlin), dwork(pdccs),
     &           dwork(pdbcs), dwork(pdinf+1), dwork(pdsup+1),
     &           neqtot, nintot,
     &           iwork(piw), dwork(pdw), dwork(pdlagr), wopt,
     &           info )
c      
      RETURN
      END
