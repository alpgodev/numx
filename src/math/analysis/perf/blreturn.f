c=======================================================================
c     subroutine BLRETURN                                      
c
c     This function computes the expected return by the Black-Litterman 
c     method
c     
c     mu = inv(inv(taux*cov) + P'*inv(omega)*P)*(inv(alpha*cov)*rho + P'*inv(omega)*q);
c
c     mu    : expected returns
c     alpha : confidence on the equilibrium returns
c     cov   : covariance matrix
c     rho   : equilibrium returns = delta*cov*w
c     pview : matrix of investor views
c     omega : uncertainty view matrix
c     view  : investor views
c 
c----------------------------------------------------------------------
      SUBROUTINE blreturn ( p, q, alpha, cov, rho, pview, omega, view,
     &                      iwork, dwork, mu, info )
c----------------------------------------------------------------------
c
c     INPUT
c       p     : portfolio size (p>1)                           integer
c       q     : number of view (q>0)                            double
c       alpha : confidence on the equilibrium returns           double
c       cov   : covariance matrix (p by p)                      double
c       rho   : equilibrium returns = delta*cov*w (p)           double
c       pview : matrix of investor views (q by p)               double
c       omega : uncertainty view matrix  (q by q)               double
c       view  : investor views (q)
c
c     WORKSPACES
c       iwork :                                                 double
c       dwork :                                                 double
c
c     OUTPUT 
c       mu    : Black-Litterman estimator (p)                   double
c       info  : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, q, info
      DOUBLE PRECISION mu(*), rho(*), cov(*), pview(*), omega(*), 
     &                 view(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, pdones, pdinvC, pdjms, pdv
      DOUBLE PRECISION x, y, t, w, EPS, ZERO, ONE
      PARAMETER ( EPS = 1.E-30, ZERO = 0.0, ONE = 1.0 )
c
c     external subroutines
      EXTERNAL PMX, JMS, OVTMCV
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      DO i = 1,p
        mu(i) = rho(i)
      ENDDO
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------     
      pdjms  = 1
c     pdjms  : pointer for workspaces of JMS (p*p)            
      pdinvC = pdjms + ( p*p )
c     pdinvC : pointer for inverse cov. matrix (p*p)
      pdones = pdinvC + ( p*p )
c     pdones : pointer for ones vector (p)
      pdv   = pdones + ( p )
c     pdv   : pointer for temporary vector                          
c
c     Total size of dwork array = 2*p*(p + 1)
c
c----------------------------------------------------------------------
c
c     alpha*cov
      CALL PMX ( p, p, cov, alpha, dwork(pdcov)) )
c
c     inv(alpha*cov)
      CALL JMS ( p, dwork(pdcov), dwork(pdjms), dwork(pdinvC), info ) 
      IF (info .NE. 0) THEN
        info = -108
        RETURN
      ENDIF
c           
c     P'*inv(omega)*P
      CALL OVTMCV (p, dwork(pdinvC), dwork(pdones), y )
      IF (y .LT. EPS) THEN
        info = -1
        RETURN
      ENDIF
c     1'*inv(cov)*rho
c     matrix vector multiply: U(p) = A(p,p)*V(p)
      CALL PMV ( p, p, dwork(pdinvC), rho, dwork(pdv) )
      CALL SEV ( p, dwork(pdv), x )
c
      RETURN
      END
