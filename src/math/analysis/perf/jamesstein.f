c=======================================================================
c     subroutine JAMESSTREIN                                      
c
c     James-Stein expected return - shrinkage estimator
c     
c     mu = (1 - intensity)*mu + intensity*target
c
c     mu        : expected return
c     intensity : shrinkage intensity
c     target    : shrinkage taget 
c 
c----------------------------------------------------------------------
      SUBROUTINE jamesstein ( n, p, cov, rho, iwork, dwork, mu, info )
c----------------------------------------------------------------------
c
c     INPUT
c       n     : number of observations (n>1)                   integer  
c       p     : portfolio size (p>1)                           integer
c       cov   : covariance matrix (p*p)                         double
c       rho   : expected return (p)                             double
c
c     WORKSPACES
c       iwork : 12*p + 3                                       integer 
c       dwork : p*(19*p + 63)/2 + 2                             double
c
c     OUTPUT 
c       mu    : James-Stein estimator (p)                       double
c       info  : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION mu(*), rho(*), cov(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
      INTEGER iwork(*)
c
c     local variables
      INTEGER i, mcte, pisdls, pdcov, pdsdls, pdA, pdb, pdinvC, pdjms, 
     &        pdv, pdones
      DOUBLE PRECISION x, y, t, w, epsbfg, alpha, EPS, ZERO, ONE, m
      PARAMETER ( EPS = 1.E-30, ZERO = 0.0, ONE = 1.0 )
c
c     external subroutines
      EXTERNAL sdls, IVX, JMS, OVTMCV, PMV, SEV
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      CALL IVX ( p, mu, ZERO )
c
c     pointers for integer work space  : iwork
c     ----------------------------------------     
      pisdls = 1
c     pisdls : pointer for SDLS workspaces (12*p + 3)                     
c
c     Total size of dwork array = 12*p + 3
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------     
      pdcov  = 1
c     pdcov  : pointer for temporary cov. matrix (p*p)
      pdsdls = pdcov + ( p*p )
c     pdsdls : pointer for SDLS workspaces (p*(p+1)*5/2 + p*(4*p + 27))
      pdA    = pdsdls + ( p*(p+1)*5/2 + p*(4*p + 27) )
c     pdA    : pointer for constraint matrix (1)
      pdb    = pdA + ( 1 )
c     pdb    : pointer for constraint vector (1)    
      pdjms  = pdb + ( 1 )
c     pdjms  : pointer for JMS workspaces (p*p)            
      pdinvC = pdjms + ( p*p )
c     pdinvC : pointer for inverse cov. matrix (p*p)
      pdones = pdinvC + ( p*p )
c     pdones : pointer for ones vector (p)
      pdv   = pdones + ( p )
c     pdv   : pointer for temporary vector (p)                         
c
c     Total size of dwork array = 2*p*(p + 1)
c                               = p*(19*p + 63)/2 + 2 
c
c----------------------------------------------------------------------
c
c     shrinkage target
c     inv(cov)
      CALL IVX ( p, dwork(pdones), ONE )
      CALL JMS ( p, cov, dwork(pdjms), dwork(pdinvC), info )
c
c     matrix correction to ensure definite positivity
      IF (info .NE. 0) THEN
        mcte = 0        ! number of constraints
        epsbfg = 1.E-8  ! BFGS precision stop test
        alpha  = 1.E-8  ! minimum eigenvalue
        CALL sdls ( p, cov, mcte, dwork(pdA), dwork(pdb), epsbfg, alpha,
     &              iwork(pisdls), dwork(pdsdls), dwork(pdcov), info )
        CALL JMS ( p, dwork(pdcov), dwork(pdjms), dwork(pdinvC), info )
      ENDIF
c           
c     1'*inv(cov)*1
      CALL OVTMCV (p, dwork(pdinvC), dwork(pdones), y )
      IF (y .LT. EPS) THEN
        info = -1
        RETURN
      ENDIF
c      
c     1'*inv(cov)*rho
c     matrix vector multiply: U(p) = A(p,p)*V(p) 
      CALL PMV ( p, p, dwork(pdinvC), rho, dwork(pdv) )
      CALL SEV ( p, dwork(pdv), x )
c
c     global minimum variance portfolio
c     target (1'*inv(cov)*rho)/(1'*inv(cov)*1)     
      t = x/y 
c
c     shrinkage intensity 
c     (p+2)/(p+2+(n*(n-1)/(n-p-2))*(rho-t)'*inv(cov)*(rho-t))
      DO i = 1,p
        dwork(pdv + i - 1) = rho(i) - t
      ENDDO
      CALL OVTMCV (p, dwork(pdinvC), dwork(pdv), x )
      IF ((N-p-2.0).LT.EPS) THEN
        info = -1
        RETURN
      ENDIF
      m = (N*(N-1.0))/(N-p-2.0)
      w = (p + 2.0)/(p + 2.0 + m*x)
c
c     James-Stein shrinkage estimator
      DO i = 1,p
        mu(i) = (1 - w)*rho(i) + w*t
      ENDDO  
      RETURN
      END
