c=======================================================================
c
c     Asset Allocation Utilities Functions                   
c
c-----------------------------------------------------------------------
c
c        CTIT    : constraints modelling (ALLOCIT)
c        CTITCOV : quadratic part (ALLOCMVTC)
c        CTITCST : constraints modelling (ALLOCMVTC)
c        CTMV    : constraints modelling (ALLOCMV)
c        CTSR    : constraints modelling (ALLOCSR)
c        CTVOL   : constraints modelling (ALLOCVOL)
c        CTMVRFR : constraints modelling (ALLOCMVRFR)
c
c=======================================================================
c
c     subroutine CTIT
c
c     Computing constraints for ALLOCIT - Index Tracking allocation
c
c-----------------------------------------------------------------------
      SUBROUTINE ctit ( n, rhob, rho, covb,
     &                  delta, neq, nin, ccst, bcst,
     &                  neqtot, nintot, ccstot, bcstot, plin )    
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                             integer
c            rhob   : mean return of index                        double
c            rho    : mean returns (n)                            double
c            covb   : covariance assets/index (n)                 double
c            delta  : outperformance target                       double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : constraints matrix (n*(neq + nin))          double
c            bcst   : constraints vector (neq + nin)              double
c
c     OUTPUT 
c            neqtot : number of equality constraints             integer
c            nintot : number of inequality constraints           integer
c            ccstot : constraints matrix (n*(neq + nin + 1))      double
c            bcstot : constraints vector (neq + nin + 1)          double
c            plin   : linear part (n)                             double
c
c     METHOD 
c                 composition of ccstot, bcstot, plin
c         - ccstot
c           | neq | nin |      1      |  size
c            _________________________
c           |           |  -rho(1)    |
c           |           |     ..      |
c           |           |             |
c           |    ccst   |             |     n
c           |           |             |  ( w(i) )
c           |           |             |
c           |           |  -rho(n)    |
c           __________________________
c
c         - bcstot
c           | neq | nin |      1      |  size
c           __________________________
c           |    bcst   | -delta-rhob |   1
c           __________________________
c
c         - plin
c           |          n            |  size
c           ________________________
c           |      -covb            |   1
c           ________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION rhob, delta
      DOUBLE PRECISION rho(*), covb(*), bcst(*), bcstot(*), plin(*),
     &                 ccst(n,*), ccstot(n,*)
c
c     local variables
      INTEGER i, j, ncst
c
c-----------------------------------------------------------------------
c
c     initializations 
      ncst = neq + nin   
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + 1
c
c     linear part
      DO i = 1,n
         plin(i) = -covb(i)
      ENDDO
c
c     case ncst = 0 (neq=0, nin=0)
      IF (ncst .EQ. 0) THEN
         DO i = 1,n
           ccstot(i, 1) = -rho(i)          
         ENDDO
         bcstot(1) = -delta -rhob
         RETURN
      ENDIF         
c
c     constraints matrix C*w <= b 
      DO i = 1,n
         DO j = 1,ncst
            ccstot(i, j) = ccst(i, j)
         ENDDO
         ccstot(i, ncst+1) = -rho(i)
      ENDDO
c
c     constraints vector
      DO j = 1,ncst
         bcstot(j) = bcst(j)
      ENDDO
      bcstot(ncst+1) = -delta - rhob
c
      RETURN
      END
c      
c=======================================================================
c
c     subroutine CTITCST
c
c     Computing constraints (portfolio allocation with transactions costs)
c
c-----------------------------------------------------------------------
      SUBROUTINE ctitcst ( n, cost, neq, nin, ccst, bcst, wini,
     &                     neqtot, nintot, ccstot, bcstot, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of asset(s)                         integer
c            cost   : transaction cost (n)                        double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix initial of constraints (n*(neq+nin)) double
c            bcst   : vector initial of constraints (neq+nin)     double
c            wini   : initial weight(s) (n)                       double
c
c     OUTPUT 
c            neqtot : number of equality constraints             integer
c            nintot : number of inequality constraints           integer
c            ccstot : matrix of constraints (2n*(neq+nin+2n))     double
c            bcstot : vector of constraints (neq+nin+2n)          double
c            info   : diagnostic argument                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot, info
      DOUBLE PRECISION bcst(*), bcstot(*), wini(*), ccst(n,*), 
     &                 ccstot(2*n,*), cost(*)
c
c     local variables
      INTEGER p, i, j, ncst, ncstot
      DOUBLE PRECISION a, zero
      PARAMETER (zero = 0.0)
c      INTEGER unit  
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
      p    = n + n
      ncst = neq + nin
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + n + n
c
c     total nb of constraints (C + cost transaction)      
      ncstot = neqtot + nintot      
c
c     constraints matrix initialization
      CALL IMX ( p, ncstot, ccstot, zero )
c      DO i = 1,p
c         DO j = 1,ncstot
c            ccstot(i , j) = 0.0
c         ENDDO
c      ENDDO
c
c     constraint C*w <= b 
      IF (ncst .NE. 0) THEN
        DO i = 1,n
            DO j = 1,ncst
                ccstot(i, j) = ccst(i, j)
            ENDDO
        ENDDO  
        DO j = 1,ncst
            bcstot(j) = bcst(j)
        ENDDO
      ENDIF  
c
c     transaction cost constraints
      DO i = 1,n
        ccstot(i, ncst + i)         =  1.0
        ccstot(i + n, ncst + i)     = -1.0
        ccstot(i, ncst + n + i)     = -1.0
        ccstot(i + n, ncst + n + i) = -1.0
        
        bcstot(ncst + i)            =  wini(i)
        bcstot(ncst + n + i)        = -wini(i)
      ENDDO
      a=cost(1)
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine CTITCOV
c
c     Computing quadratic part for ALLOCMVTC
c     Index Tracking allocation with cost transaction
c
c-----------------------------------------------------------------------
      SUBROUTINE ctitcov ( n, cov, alpha, cost, gamma, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of asset(s)                         integer
c            cov    : covariance matrix (n*n)                     double
c            alpha  : cost transaction parameter                  double
c            cost   : cost transaction (n)                        double
c
c     OUTPUT 
c            gamma  : new quadratic part (2n*2n)                  double
c            info   : diagnostic argument                        integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha
      DOUBLE PRECISION cost(*), cov(n,*), gamma(n + n,*)
c
c     local variables
      INTEGER p, i, j
      DOUBLE PRECISION zero, eps
      PARAMETER (zero = 0.0, eps = 1.E-5)
c      INTEGER ioecr 
c
c-----------------------------------------------------------------------
c
c     initialisations
      info = 0 
      p    = n + n
c
c     construction of the quadratic part (matrix)
c
c     gamma = [cov        0        ] 
c             [ 0   (alpha*cost)^2 ]
      CALL IMX ( p, p, gamma, zero )
c      DO i = 1,p
c        DO j = 1,p
c            gamma(i, j) = 0.0
c        ENDDO
c      ENDDO
      DO i = 1,n
        DO j = 1,n
            gamma(i, j) = cov(i, j)
        ENDDO
        gamma(n+i, n+i) = max(alpha*alpha*cost(i)*cost(i), n*eps)
      ENDDO
c
      RETURN
      END
c
c=======================================================================      
c
c     subroutine CTMV
c
c     Computing constraints for ALLOCMV - Mean Variance (Markowitz)
c
c-----------------------------------------------------------------------
      SUBROUTINE ctmv ( n, rho, mu, neq, nin, ccst, bcst,
     &                  neqtot, nintot, ccstot, bcstot, plin )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                             integer
c            rho    : mean returns (n)                            double
c            mu     : performance target                          double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : constraints matrix (n*(neq + nin))          double
c            bcst   : constraints vector (neq + nin)              double
c
c     OUTPUT 
c            neqtot : number of equality constraints             integer
c            nintot : number of inequality constraints           integer
c            ccstot : matrix of constraints (n*(neq+nin+1))       double
c            bcstot : vector of constraints (neq+nin+1)           double
c            plin   : linear part (n)                             double
c
c     METHOD :
c                 composition of ccstot, bcstot, plin
c         - ccstot
c           |     neq   |   nin   | 1 |  size
c            _________________________
c           |                     |   |
c           |                     |   |
c           |                     | - |
c           |          ccst       | r |     n
c           |                     | h |  ( w(i) )
c           |                     | o |
c           |                     |   |
c           __________________________
c
c         - bcstot
c           |     neq   |   nin   | 1 |  size
c           __________________________
c           |          bcst       |-mu|   1
c           __________________________
c
c         - plin
c           |          n            |  size
c           ________________________
c           |            0          |   1
c           ________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION mu
      DOUBLE PRECISION rho(*), bcst(*), bcstot(*), plin(*), 
     &                 ccst(n,*), ccstot(n,*)
c
c     local variables
      INTEGER i, j, ncst
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.D0)
c
c-----------------------------------------------------------------------
c
c     initialisations 
      ncst = neq + nin
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + 1
c
c     construction of the linear part
      CALL IVX ( n, plin, ZERO )
c      
c     case ncst = 0 (neq=0, nin=0)
      IF (ncst .EQ. 0) THEN
         DO i = 1,n
            ccstot(i, 1) = -rho(i)
         ENDDO
         bcstot(1) = -mu
         RETURN
      ENDIF
c
c     constraints matrix
      DO i = 1,n
         DO j = 1,ncst
            ccstot(i,j) = ccst(i,j)
         ENDDO
         ccstot(i,ncst+1) = -rho(i)
      ENDDO
c
c     constraints vector
      DO j = 1,ncst
        bcstot(j) = bcst(j)
      ENDDO  
      bcstot(ncst+1) = -mu
c
      RETURN
      END
c
c=======================================================================     
c
c     subroutine CTSR
c
c     Utility for OPSR (Max. Sharpe Ratio)
c
c-----------------------------------------------------------------------
      SUBROUTINE ctsr ( n, rho, rfr, neq, nin, ccst, bcst, cinf, csup,
     &                  neqtot, nintot, ccstot, bcstot)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size of portfolio                               integer
c       rho    : mean returns (n)                                 double
c       rfr    : risk-free rate                                   double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : matrix of constraints (n*(neq+nin))              double
c       bcst   : vector of constraints (neq+nin)                  double
c       cinf   : lower bounds (n)                                 double
c       csup   : upper bounds (n)                                 double
c
c     OUTPUT 
c       neqtot : total number of equality constraints            integer
c       nintot : total number of inequality constraints          integer
c       ccstot : matrix of constraints 
c               (n+1)*(2+2*n+nin)                                 double
c       bcstot : vector of constraints (2+2*n+nin)                double
c
c     METHOD 
c
c         - ccstot
c           |(1)|   (1)  | (n)  |  (n)   |(nin) |    (size)
c            ____________________________________
c           | 1 | rho(1) |      |        |      |   
c           | . | .      | -Id  |   Id   |  cst |   (n + 1)
c           | 1 | rho(n) |      |        |      |
c           |------------------------------------
c           |-1 | -rfr   | Cinf | -Csup  |  -b  |     (1)
c            ____________________________________
c
c         - bcstot
c           |(1)|   (1)  | (n)  |  (n)   |(nin) |    (size)
c           _____________________________________ 
c           | 0 |    1   |  0   |   0    |  0   |     (1)
c           _____________________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION rfr
      DOUBLE PRECISION rho(*), bcst(*), bcstot(*), cinf(*), csup(*), 
     &                 ccst(n,*), ccstot(n+1,*)
c
c     local variables
      INTEGER i, j, p, ncst
c
c-----------------------------------------------------------------------
c
c     number of variables
      p = n + 1
c
c     nb of equality and inequality constraints
      neqtot = 2
      nintot = nin + 2*n
c
c     total number of constaints      
      ncst   = neqtot + nintot
c
c     construction of the constraints matrix
      DO i = 1,n
         ccstot(i, 1) = 1
         ccstot(i, 2) = rho(i)
         IF (nin .GT. 0) THEN
            DO j = 1,nin
                ccstot(i, neqtot + j) = ccst(i,neq + j)
            ENDDO
         ENDIF
         DO j = 1,n
            ccstot(i, neqtot + nin + j)     = 0.
            ccstot(i, neqtot + nin + n + j) = 0.
            IF (i .EQ. j) THEN
                ccstot(i, neqtot + nin + j)     = -1.0
                ccstot(i, neqtot + nin + n + j) =  1.0
            ENDIF
         ENDDO
      ENDDO
c
      ccstot(p, 1) = -1.0
      ccstot(p, 2) = -rfr
      IF (nin .GT. 0) THEN
        DO j = 1,nin
            ccstot(p, neqtot + j) = -bcst(neq + j)
        ENDDO
      ENDIF
      DO j = 1,n
        ccstot(p, neqtot + nin + j)     =  cinf(j)
        ccstot(p, neqtot + nin + n + j) = -csup(j)
      ENDDO
c
c     construction of the constraints vector
      bcstot(1) = 0.0
      bcstot(2) = 1.0
      IF (nin .GT. 0) THEN
        DO j = 1,nin
            bcstot(neqtot + j) = 0.
        ENDDO
      ENDIF
      DO j = 1,(2*n)
        bcstot(neqtot + nin + j) = 0.
      ENDDO
c
      RETURN
      END
c
c
c=======================================================================      
c
c     subroutine CTVOL
c
c     Computing constraints for ALLOCVOL
c
c-----------------------------------------------------------------------
      SUBROUTINE ctvol ( n, rho, mu, neq, nin, ccst, bcst,
     &                   neqtot, nintot, ccstot, bcstot, plin )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : portfolio size                             integer
c            rho    : mean returns (n)                            double
c            mu     : performance target                          double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : constraints matrix (n*(neq + nin))          double
c            bcst   : constraints vector (neq + nin)              double
c
c     OUTPUT 
c            neqtot : number of equality constraints             integer
c            nintot : number of inequality constraints           integer
c            ccstot : matrix of constraints (n*(neq+nin+1))       double
c            bcstot : vector of constraints (neq+nin+1)           double
c            plin   : linear part (n)                             double
c
c     METHOD :
c                 composition of ccstot, bcstot, plin
c         - ccstot
c           | 1 |    neq   |   nin   | 1 |  size
c            __________________________
c           | 1 |                  |   |
c           |   |                  |   |
c           | . |                  | - |
c           | . |       ccst       | r |     n
c           | . |                  | h |  ( w(i) )
c           |   |                  | o |
c           | 1 |                  |   |
c           ___________________________
c
c         - bcstot
c           | 1 |    neq   |   nin   | 1 |  size
c           _____________________________
c           | 1 |         bcst       |-mu|   1
c           _____________________________
c
c         - plin
c           |          n            |  size
c           ________________________
c           |          0            |   1
c           ________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION mu
      DOUBLE PRECISION rho(*), bcst(*), bcstot(*), plin(*), 
     &                 ccst(n,*), ccstot(n,*)
c
c     local variables
      INTEGER i, j, ncst, ncsttot
c
c-----------------------------------------------------------------------
c
c     initialisations 
      ncst = neq + nin
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + 1
      ncsttot = neqtot + nintot
c
c     construction of the linear part (equal zero)
      DO i = 1,n
         plin(i) = 0.0
      ENDDO
c
c     constraints: w'*rho >= mu            
      IF (ncst .EQ. 0) THEN
c        case ncst = 0 (neq=0, nin=0)
         DO i = 1,n
            ccstot(i, 1) = -rho(i)
         ENDDO
         bcstot(1) = -mu
         RETURN
      ENDIF   
c
c     general case ncst <> 0
      DO i = 1,n
        DO j = 1,ncst
            ccstot(i, j) = ccst(i, j)
        ENDDO
        ccstot(i,ncst+1) = -rho(i)
      ENDDO
c
c     constraints vector
      DO j = 1,ncst
        bcstot(j) = bcst(j)
      ENDDO  
      bcstot(ncst+1) = -mu     
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine CTMVRFR
c
c     Computing constraints with risk-free rate
c
c-----------------------------------------------------------------------
      SUBROUTINE ctmvrfr ( n, rho, rfr, mu, neq, nin, ccst, bcst,
     &                     cinfrfr, csuprfr,
     &                     neqtot, nintot, ccstot, bcstot, plin )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       rho    : expected mean returns (n)                        double
c       rfr    : risk-free rate                                   double
c       mu     : performance target                               double
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       ccst   : constraints matrix (n*(neq + nin))               double
c       bcst   : constraints vector (neq + nin)                   double
c       cinfrfr: risk-free rate lower bound                       double
c       csuprfr: risk-free rate upper bound                       double
c
c     OUTPUT 
c       neqtot : number of equality constraints                  integer
c       nintot : number of inequality constraints                integer
c       ccstot : matrix of constraints (n*(neq+nin+3))            double
c       bcstot : vector of constraints (neq+nin+3)                double
c       plin   : linear part (n)                                  double
c
c     METHOD : composition of ccstot, bcstot, plin
c         - ccstot
c           |     neq   |   nin   | 1 | 1 |  1 | size
c            __________________________________
c           |                     |   | 1 | -1 |
c           |                     |   | . |  . |
c           |                     | - | . |  . |
c           |          ccst       | r | . |  . |  n
c           |                     | h | . |  . | 
c           |                     | o | . |  . |
c           |                     |   | 1 | -1 |
c           ____________________________________
c
c         - bcstot
c           | neq | nin | 1 |     1     |     1     |   size
c           _________________________________________
c           |   bcst    |-mu| 1-Cinfrfr | Csuprfr-1 |    1
c           _________________________________________
c
c         - plin
c           |          n            |  size
c           ________________________
c           |          0            |   1
c           ________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION mu, rfr, cinfrfr, csuprfr
      DOUBLE PRECISION rho(*), bcst(*), bcstot(*), plin(*), 
     &                 ccst(n,*), ccstot(n,*)
c
c     local variables
      INTEGER i, j, ncst
      DOUBLE PRECISION ZERO, EPS
      PARAMETER (ZERO = 0.D0, EPS = 1.E-8)
c
c-----------------------------------------------------------------------
c
c     initialisations 
      ncst = neq + nin
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + 3
c
c     construction of the linear part
      CALL IVX (n, plin, ZERO)

c     case ncst <> 0 (neq=0, nin=0)
c      IF (ncst .NE. 0) THEN
c        DO i = 1,n
c            DO j = 1,ncst
c                ccstot(i, j) = ccst(i, j)
c            ENDDO
c        ENDDO  
c        DO j = 1,ncst
c            bcstot(j) = bcst(j)
c        ENDDO
c      ENDIF
c     matrix constraint      
c      DO i = 1,n
c        ccstot(i, ncst + 1) =  rfr - rho(i)
c        ccstot(i, ncst + 2) =  1.0
c        ccstot(i, ncst + 3) = -1.0
c      ENDDO
c        bcstot(ncst + 1) = rfr - mu
c        bcstot(ncst + 2) = 1.0 - cinfrfr + EPS
c        bcstot(ncst + 3) = csuprfr - 1.0 + EPS
c
c      
c     case ncst = 0 (neq=0, nin=0)
      IF (ncst .EQ. 0) THEN
         DO i = 1,n
            ccstot(i, 1) =  rfr - rho(i)
            ccstot(i, 2) =  1.0
            ccstot(i, 3) = -1.0
         ENDDO
         bcstot(1) = rfr - mu
         bcstot(2) = 1.0 - cinfrfr + EPS
         bcstot(3) = csuprfr - 1.0 + EPS
         RETURN
      ENDIF
c
c     constraints matrix
      DO i = 1,n
         DO j = 1,ncst
            ccstot(i, j) = ccst(i, j)
         ENDDO
         ccstot(i,ncst + 1) =  rfr - rho(i)
         ccstot(i,ncst + 2) =  1.0
         ccstot(i,ncst + 3) = -1.0
      ENDDO
c
c     constraints vector
      DO j = 1,ncst
        bcstot(j) = bcst(j)
      ENDDO  
      bcstot(ncst+1) = rfr - mu
      bcstot(ncst+2) = 1.0 - cinfrfr + EPS
      bcstot(ncst+3) = csuprfr - 1.0 + EPS
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine CTITRFR
c
c     Computing constraints for ALLOCITRFR - Index Tracking allocation
c
c-----------------------------------------------------------------------
      SUBROUTINE ctitrfr ( n, rhob, rho, covb,
     &                     delta, neq, nin, ccst, bcst, 
     &                     cinfrfr, csuprfr,
     &                     neqtot, nintot, ccstot, bcstot, plin )    
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : portfolio size                                 integer
c       rhob    : mean return of index                            double
c       rho     : mean returns (n)                                double
c       covb    : covariance assets/index (n)                     double
c       delta   : outperformance target                           double
c       neq     : number of initial equality constraints         integer
c       nin     : number of initial inequality constraints       integer
c       ccst    : constraints matrix (n*(neq + nin))              double
c       bcst    : constraints vector (neq + nin)                  double
c       cinfrfr : risk-free rate lower bound                      double
c       csuprfr : risk-free rate upper bound                      double
c
c     OUTPUT 
c       neqtot  : number of equality constraints                 integer
c       nintot  : number of inequality constraints               integer
c       ccstot  : constraints matrix (n*(neq + nin + 3))          double
c       bcstot  : constraints vector (neq + nin + 3)              double
c       plin    : linear part (n)                                 double
c
c     METHOD 
c                 composition of ccstot, bcstot, plin
c         - ccstot
c           | neq | nin |      1      |  size
c            _________________________
c           |           |  -rho(1)    |
c           |           |     ..      |
c           |           |             |
c           |    ccst   |             |     n
c           |           |             |  ( w(i) )
c           |           |             |
c           |           |  -rho(n)    |
c           __________________________
c
c         - bcstot
c           | neq | nin |      1      |  size
c           __________________________
c           |    bcst   | -delta-rhob |   1
c           __________________________
c
c         - plin
c           |          n            |  size
c           ________________________
c           |      -covb            |   1
c           ________________________
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, neqtot, nintot
      DOUBLE PRECISION rhob, delta, cinfrfr, csuprfr
      DOUBLE PRECISION rho(*), covb(*), bcst(*), bcstot(*), plin(*),
     &                 ccst(n,*), ccstot(n,*)
c
c     local variables
      INTEGER i, j, ncst
      DOUBLE PRECISION ZERO, EPS
      PARAMETER (ZERO = 0.D0, EPS = 1.E-8)
c
c-----------------------------------------------------------------------
c
c     initializations 
      ncst = neq + nin
c
c     nb of equality and inequality constraints
      neqtot = neq
      nintot = nin + 3
c
c     linear part
      DO i = 1,n
         plin(i) = -covb(i)
      ENDDO
c      
c     case ncst <> 0 (neq=0, nin=0)
      IF (ncst .NE. 0) THEN
        DO i = 1,n
            DO j = 1,ncst
                ccstot(i, j) = ccst(i, j)
            ENDDO
        ENDDO  
        DO j = 1,ncst
            bcstot(j) = bcst(j)
        ENDDO
c        CALL YM ( n, ncst, ccst, ccstot ) ! constraint matrix
c        CALL YV ( ncst, bcst, bcstot )    ! constraint vector
      ENDIF
c     matrix constraint      
      DO i = 1,n
        ccstot(i, ncst + 1)=  -rho(i)
        ccstot(i, ncst + 2) =  1.0
        ccstot(i, ncst + 3) = -1.0
      ENDDO
      bcstot(ncst + 1) = -delta - rhob
      bcstot(ncst + 2) = 1.0 - cinfrfr + EPS
      bcstot(ncst + 3) = csuprfr - 1.0 + EPS
c
      RETURN
      END
