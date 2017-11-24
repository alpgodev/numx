c=======================================================================
c
c     subroutine APTCST                                     
c
c     APT (Arbitrage Pricing Theory) model with alpha/beta constraints
c
c            Y = alpha + X*beta'
c 
c     where Y is the matrix (n*p) of asset(s) value(s), X the matrix (n*q) 
c     of factor(s) value(s) and beta=[beta(1), ...,beta(q)]
c     
c     Min[ theta'*(F'F)*theta - 2*theta'*(F'*Y) ], 
c     where F = [1., X] and theta = [alpha, beta']
c      s.t. 
c     C*theta <= b             (linear constraints) 
c     Cinf <= theta <= Csup    (lower/upper bounds)
c
c-----------------------------------------------------------------------
      SUBROUTINE aptcst ( n, q, Y, X,
     &                    neq, nin, ccst, bcst, cinf, csup,
     &                    iwork, dwork,
     &                    alpha, beta, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of value(s) (n > 0)                 integer
c            q      : number of factors  (q > 0)                 integer
c            Y      : asset(s) values (n)                         double
c            X      : factor(s) values (n*q)                      double
c            neq    : number equality constraints                integer
c            nin    : number inequality constraints              integer
c            ccst   : matrix of constraints ((q+1)*(neq+nin))     double
c            bcst   : vector initial of constraints (neq+nin)     double
c            cinf   : lower bound (q+1)                           double
c            csup   : upper bound (q+1)                           double
c
c     WORKSPACE 
c            iwork  : 3*p + 2*nin + neq + 1                      integer 
c            dwork  : p*(2*p + n + 9) + 3*nin + neq               double
c                     with p = q + 1 
c
c     OUTPUT 
c            alpha  : alpha coefficient(s)                        double
c            beta   : beta coefficient(s) (q)                     double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            IVX    : initialization at a scalar of a vector
c            YV     : copy a vector in a vector
c            PRMM   : M'*M -> Matrix
c            PMTV   : M'*V -> Vector
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, q, info, neq, nin
      DOUBLE PRECISION alpha
      DOUBLE PRECISION Y(*), X(*), cinf(*), csup(*), ccst(q,*), bcst(*),
     &                 beta(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, p, piw, pf, pdqua, pdlin, pdw, pdlagr, pdxopt
      DOUBLE PRECISION DUN
      PARAMETER ( DUN = 1.E0 )
c
c     external subroutines
      EXTERNAL PRMM, PMTV, YV, IVX, qp
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
      p = q + 1
c
c     pointers for integer work space : iwork
c     ---------------------------------------  
      piw  = 1
c     piw  : pointer for internal workspaces of QUAPRO
c             needs( 3*p + 2*nin + neq + 1 )     
c
c     Total size of iwork array = 3*p + 2*nin + neq + 1
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pf     = 1
c     pf    : pointer for matrix F, so ( n*p ) more 
      pdqua  = pf + ( n*p )
c     pdqua : pointer for quadratic part matrix( p*p ) 
      pdlin  = pdqua + ( p*p )
c     pdlin : pointer for linear part, so ( p ) more
      pdw    = pdlin + ( p )
c     pdw   : pointer for QUAPRO, so (p*p + 6*p + 2*nin) more
      pdlagr = pdw + ( p*p + 6*p + 2*nin )
c     pdlagr: pointer for QUAPRO Lagrange multipliers, so (p+nin+neq) more
      pdxopt = pdlagr + ( p + nin + neq )
c     pdxopt: pointer for QUAPRO optimal points, ( p )
c
c     Total size of dwork array = n*p
c                               + p*p 
c                               + p 
c                               + p*p + 6*p + 2*nin 
c                               + p + nin + neq
c                               + p
c                              =  p*(2*p + n + 9) + 3*nin + neq
c
c-----------------------------------------------------------------------
c
c     construction of matrix F = [1.0, X]
      CALL IVX ( n, dwork(pf), DUN )
      CALL YV ( q*n, X, dwork(pf + n) )
c
c     quadratic part, F'*F
      CALL PRMM ( n, p, dwork(pf), dwork(pdqua) )
c
c     test matrix condition number (Lmax/Lmin)
      
c
c     applied SDLS correction if bad condition number
      
c
c     linear part, -2*(F'*Y) 
      CALL PMTV ( n, p, dwork(pf), Y, dwork(pdlin) )
      DO i = 1,p
        dwork(pdlin + i - 1) = -2.0*dwork(pdlin + i - 1)
      ENDDO
c
c     quadratic solver: Min[ theta'*(F'F)*theta - 2*theta'*(F'*Y) ]
      CALL qp ( p, dwork(pdqua), dwork(pdlin), neq, nin, ccst,
     &          bcst, cinf, csup,
     &          iwork(piw), dwork(pdw), dwork(pdlagr),dwork(pdxopt),
     &          info )
c
c     coefficients alpha and beta
      alpha = dwork(pdxopt)
      DO i = 1,q
        beta(i) = dwork(pdxopt + i)
      ENDDO
c
      RETURN
      END
