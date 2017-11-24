c=======================================================================
c
c     subroutine BUILDRBIT                                     
c
c     Uses the inputs of riskbudgetIT to compute the suitable inputs for
c     the solver riskbudegting
c
c-----------------------------------------------------------------------
      SUBROUTINE buildrbit ( n, cov, kappa, rho, covb, varb,
     &                       neq, nineq, ccst, bcst, cinf, csup,
     &                       nbc, class, sigma,
     &                       iwork, dwork, omegabench,
     &                       rhoIT, bcstIT, cinfIT, csupIT, sigmaIT,
     &                       info )
c-----------------------------------------------------------------------
c
c     INPUT 
c
c       n      : portfolio size                                  integer
c       cov    : covariance matrix (n*n)                          double
c       rho    : expected returns vector (n)                      double
c       covb   : covariance assets-index (n)                      double
c       varb   : variance of the benchmark                        double
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       ccst   : matrix of constraints (nasset*(neq+nin))         double
c       bcst   : vector initial of constraints (neq+nin)          double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       nbc    : number of risk budgeting constraints (>=1)      integer 
c       class  : block definition (n)                            integer
c       sigma  : volatilities budgets constraints (nbc)           double 
c
c     WORKSPACE 
c       iwork  : 12*n + nbc                                        integer 
c       dwork  : n*(4*n+28) + nbc + nineq + neq                     double
c                    
c     OUTPUT 
c       omegabench:  variable omega benchmark (n)                   double
c       rhoIT: n-double, equal to -(rho + kappa* Gamma_b)           double     
c       bcstIT: (neq-nineq)-double 
c       binfIT   : modified lower bound (n)                         double
c       bsupIT   : modifief upper bound (n)                         double
c       sigmaIT  : nbc-vector of Index tracking risk constraints    double
c
c     CALL   
c        IVI, IMX, SCHURS, OMCDMCT, PMV, OMEGAB    
c
c-----------------------------------------------------------------------   
c
      IMPLICIT NONE 
c
c     arguments i/o
      INTEGER n, info, neq, nineq, nbc, class(*)
      DOUBLE PRECISION cov(*), covb(*), rhoIT(*), bcstIT(*), cinfIT(*)
      DOUBLE PRECISION sigmaIT(*), csupIT(*), sigma(*), varb
      DOUBLE PRECISION kappa, rho(*), ccst(*), bcst(*), cinf(*), csup(*)
      DOUBLE PRECISION omegabench(*)
c      
c     workspaces      
      DOUBLE PRECISION dwork(*) 
      INTEGER iwork(*)
c
c     local variables
      INTEGER i, piw, pdsigmabench, pdw, pdbcstIT
      DOUBLE PRECISION NEGONE, EPS
      PARAMETER (EPS = 1.E-10, NEGONE=-1.E0)
c
c     external function
      EXTERNAL IVI, IMX, SCHURS, OMCDMCT, PMV, omegab
c
c-----------------------------------------------------------------------     
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c
c     Total size of iwork array = 12*n+nbc   
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------      
      pdsigmabench = 1
c     needs nbc
      pdbcstIT = pdsigmabench + nbc
c     needs nineq + neq
      pdw = pdbcstIT +  nineq + neq
c     needs n*(4*n+28) 
c
c     Total size of dwork array = n*(4*n+28) + nbc + nineq + neq
c      
c------------------------------------------------------------------------
c      
      CALL omegab ( n, cov, covb, nbc, class,
     &              iwork(piw), dwork(pdw),
     &              omegabench, dwork(pdsigmabench), info)
      IF (info .LT. 0) RETURN
      DO i = 1,nbc
        sigmaIT(i)=sigma(i)**2-varb+dwork(pdsigmabench+i-1)
        IF (sigmaIT(i) .LT. EPS) THEN 
            info = -112
            RETURN
        ENDIF
        sigmaIT(i) = SQRT(sigmaIT(i))
      ENDDO
      CALL PMV(n, n, cov, omegabench, rhoIT)
c
c     rhoIt = -rho + kappa( Covb + Cov*omegabench)
      IF (kappa .LT. EPS) THEN
        CALL YV(n, rho, rhoIT)
      ELSE
        DO i = 1,n
            rhoIT(i) = rhoIT(i)+covb(i)
            rhoIT(i) = -kappa*rhoIT(i)
            rhoIT(i) = rho(i)+rhoIT(i)
        ENDDO
      ENDIF
c
c     bcstIT = bcst - ccst*omegaB
      CALL PMTV(n, (neq+nineq), ccst, omegabench, dwork(pdbcstIT))
      CALL SVVX((neq+nineq), bcst, dwork(pdbcstIT), NEGONE, bcstIT)
      
c     bsupIT = bsup - omegab
      CALL SVVX(n, csup, omegabench, NEGONE, csupIT)
    
c     binfIT = binf - omegab
      CALL SVVX(n, cinf, omegabench, NEGONE, cinfIT)
      RETURN 
      END
