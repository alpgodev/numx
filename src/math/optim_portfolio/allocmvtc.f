c=======================================================================
c
c     subroutine ALLOCMVTC                                   
c
c     MV allocation srtategy with transaction costs (linear)
c
c     Min[ w*Q*w + alpha*CostTransaction(w)]
c      s.t. 
c     C*w <= b                 (linear constraints) 
c     Cinf <= w <= Csup        (lower/upper bounds)
c
c        w   : portfolio weights 
c        Q   : covariance matrix 
c        rho : assets performance 
c
c-----------------------------------------------------------------------
      SUBROUTINE allocmvtc ( n, cov, neq, nin, ccst, bcst, cinf, csup,
     &                       cost, wini, alpha,
     &                       iwork, dwork, wopt, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of assets                                integer
c       cov    : covariance matrix (n*n)                          double
c       neq    : number of equality constraints                  integer
c       nin    : number of inequality constraints                integer
c       ccst   : matrix of constraints (n*(neq+nin))              double
c       bcst   : vector of constraints (neq+nin)                  double
c       cinf   : lower bound (n)                                  double
c       csup   : upper bound (n)                                  double
c       cost   : cost transaction (n)                             double
c       wini   : initial portfolio (n)                            double
c       alpha  : transaction cost parameter                       double
c
c     WORKSPACE 
c       iwork  : 26*n + 2*nin + neq + 1                          integer 
c       dwork  : 2*n*(16*n + nin + neq + 44) + 2(2*nin + neq)     double
c 
c     OUTPUT 
c       wopt   : optimal portfolio (n)                            double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       YM      : matrix copy
c       YV      :
c       QP      : quadratic solver
c       CTITCOV :
c       CTITCST :
c       TESTSDP :
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, neq, nin
      DOUBLE PRECISION alpha
      DOUBLE PRECISION cov(*), cinf(*), csup(*), ccst(n,*), bcst(*), 
     &                 wopt(*), cost(*), wini(*)
c
c     workspaces 
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER p, i, neqtot, nintot, piwork, pdwork, piwk, 
     &        pigk, pdsp, pdss, pdsa, pdwk, piwq, pdlin, pdwq, pdlagr, 
     &        pdbcs, pdccs, pdro, pdwopt, pdcov, pdcinf, pdcsup, test
      DOUBLE PRECISION eps, epsilon, maxv
      PARAMETER ( eps = 1.E-8, epsilon = 1.E-15, maxv = 1.E+15 )
c
c     external subroutines
      EXTERNAL qp, YM, YV, ctitcov, ctitcst, testsdp
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      test = 0
      p = n
c
c     test if cost(i) = 0 for i = 1,...,n
      DO i = 1,n 
        IF (cost(i) .GT. epsilon) test = 1 
      ENDDO         
c
c     case alpha > eps (with trans. costs)      
      IF ((alpha .GT. eps).AND.(test .EQ. 1)) p = n + n
c      
c     pointers for integer workspace : iwork
c     --------------------------------------
      piwork = 1
      piwk   = piwork
c     piwk  : pointer for workspace of KATAVEP who needs ( 12*p )
      pigk   = piwk + ( 12*p )
c     pigk  : pointer Kato groups, so (p) more
c     the call of KATAVEP needs ( 13*p )
c
c     pointers for QP who uses a part of the same space
      piwq   = piwork
c     piwq  : pointer for QP who needs
c             ( 3*p + 2*nintot + neqtot + 1 ),
c          so ( 3*p + 2*(nin + p) + neq + 1 )
c          so ( 5*p + 2*nin + neq + 1 )
c
c     together KATAVEP and QP need (union of space)
c     ( 13*p + 2*nin + neq + 1 )
c
c     Total size of iwork array = ( 13*p + 2*nin + neq + 1 )
c
c     pointers for double precision workspace  : dwork
c     ------------------------------------------------
      pdwork = 1
      pdwopt  = pdwork
c     pdwopt : pointer for the optimal point, so ( p ) more    
      pdcov  = pdwopt + p
c     pdcov  : pointer for quadratic part, so (p*p) more     
      pdsp   = pdcov + ( p*p )
c     pdsp   : pointer for vector spectr of KATAVEP, so ( p ) more
      pdss   = pdsp + ( p )
c     pdss   : pointer for vector eigsor of KATAVEP, so ( p ) more
      pdsa   = pdss + ( p )
c     pdsa   : pointer for vector spectr of KATAVEP, so ( p ) more
      pdwk   = pdsa + ( p )
c     pdwk  : pointer for KATAVEP workspace, so p*(4*p + 26) more
      pdro   = pdwk + ( p*(4*p + 26) )
c     pdro   : pointer for the robust quadratic part, so ( p*p ) more  
      pdlin  = pdro + ( p*p )
c     pdlin  : pointer for linear part vector, so ( p ) more
      pdwq   = pdlin + ( p )
c     pdwq  : pointer for QP workspace who needs
c                ( p*p + 6*p + 2*nintot)
c             so ( p*p + 6*p + 2*(nin + p) )
c             so ( p*p + 8*p + 2*nin )
      pdlagr  = pdwq + ( p*p + 6*p + 2*nin )
      IF ((alpha .GT. eps).AND.(test .EQ. 1)) 
     &  pdlagr  = pdwq + ( p*p + 8*p + 2*nin )
c     pdlagr : pointer for Lagrange multipliers vector
c              ( p + nintot + neqtot) 
      pdbcs   = pdlagr + ( p + nin + neq )
      IF ((alpha .GT. eps).AND.(test .EQ. 1))
     &   pdbcs   = pdlagr + ( 2*p + nin + neq )
c     pdbcs : pointer for the constraints vector ( nintot + neqtot ),
      pdccs   = pdbcs +  ( nin + neq)
      IF ((alpha .GT. eps).AND.(test .EQ. 1))
     &   pdccs   = pdbcs +  ( p + nin + neq)
c     pdccs : pointer for the constraints matrix
c              p*( nintot + neqtot ),
      pdcinf  = pdccs + ( p*(nin + neq) )
      IF ((alpha .GT. eps).AND.(test .EQ. 1))
     &   pdcinf  = pdccs + ( p*(p + nin + neq) )
c     pdcinf : pointer for the cinf optimal point, so ( p ) more 
      pdcsup  = pdcinf + ( p )
c     pdcsup : pointer for the csup optimal point, so ( p ) more    
c
c                      = p
c                      + p*p
c                      + p
c                      + p
c                      + p
c                      + P*(4*p + 26)
c                      + p*p
c                      + p
c                      + p*( p + 8 ) + 2*nin
c                      + 2*p + nin + neq
c                      + p + nin + neq
c                      + p*( p + nin + neq )
c                      + p
c                      + p
c                      = p*( 8*p + nin + neq + 44 ) + 4*nin+ + 2*neq
c
c     Total size of dwork array = p*(8*p + nin + neq + 44) + 4*nin + 2*neq
c
c--------------------------------------------------------------------
c
c     linear part and optimal weights initialization
      DO i = 1,n
        dwork(pdlin + i - 1) = epsilon
      ENDDO   
c
c     case alpha = 0 (without transaction cost. p=0)
c     quadratic solver (Mean-Variance Optimization) 
      IF ((alpha .LE. eps).OR.(test .EQ. 0)) THEN
c
c           copy cov. matrix 
            CALL YM ( n, n, cov, dwork(pdro) )
c
c           quadratic term (matrix) SDP test
            CALL testsdp (n, dwork(pdro), iwork(piwk), dwork(pdwk),info)
            IF (info .NE. 0) THEN
                info = -108
                RETURN
            ENDIF
            CALL qp ( n, dwork(pdro), dwork(pdlin),
     &                neq, nin, ccst, bcst, cinf, csup,
     &                iwork(piwq), dwork(pdwq), dwork(pdlagr),
     &                wopt, info )
            RETURN
      ENDIF
c
c     quadratic term (matrix) SDP test
      CALL testsdp (n, cov, iwork(piwk), dwork(pdwk),info)
      IF (info .NE. 0) THEN
        info = -108
        RETURN
      ENDIF
c      
c     add transaction cost. in linear part 
c     p=[0, ..., 0, alpha*cost(1), ..., alpha*cost(n)]
      DO i = 1,n
        IF ((alpha*cost(i)) .LT. eps) THEN
            dwork(pdlin + n + i - 1) = alpha*eps
        ELSE    
            dwork(pdlin + n + i - 1) = alpha*cost(i)
        ENDIF    
      ENDDO
c
c     problem formulation: min[x'*Q*x + p'*x]
c
c     quadratic part Gamma := [cov   0             ] (cf. utalloc.f)
c                             [ 0   (alpha*cost)^2 ]
      CALL ctitcov (n, cov, alpha, cost, dwork(pdcov), info)
c
      CALL YM ( p, p, dwork(pdcov), dwork(pdro) )
c
c     quadratic term (matrix) SDP test
      CALL testsdp (p, dwork(pdro), iwork(piwk), dwork(pdwk), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF
c
c     constraint modeling (cf. UTALLOC.F)
c     C*w <= b + cost transaction
      CALL ctitcst ( n, cost, neq, nin, ccst, bcst, wini,
     &               neqtot, nintot, dwork(pdccs), dwork(pdbcs), info )
c
c     lower/upper bounds        
      DO i = 1,n
        dwork(pdcinf + i - 1)     = cinf(i)
        dwork(pdcinf + n + i - 1) = epsilon
        dwork(pdcsup + i - 1)     = csup(i)
        dwork(pdcsup + n + i - 1) = maxv
      ENDDO
c
c     quadratic solver (QUAPRO)
      CALL qp ( p, dwork(pdro), dwork(pdlin), neqtot, nintot,
     &          dwork(pdccs), dwork(pdbcs), dwork(pdcinf),
     &          dwork(pdcsup),
     &          iwork(piwq), dwork(pdwq), dwork(pdlagr),
     &          dwork(pdwopt), info )
      IF (info .LT. 0) RETURN
c
c     optimal portfolio
      CALL YV (n, dwork(pdwopt), wopt)
      RETURN
      END
