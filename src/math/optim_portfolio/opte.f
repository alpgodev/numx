c=======================================================================
c
c     subroutine OPTE
c
c     Optimization for ALLOCTE - Tracking Error Constraint
c
c-----------------------------------------------------------------------
      SUBROUTINE opte ( n, cov, rho, covb, varb, rhob,
     &                  neq, nin, ccst, bcst, cinf, csup,
     &                  temax, maxiter, epsilon,
     &                  iwork, dwork, wopt, delopt, optte, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        n      : number of asset(s) (> 1)                       integer
c        cov    : n*n-array, covariance matrix (n*n)              double  
c        rho    : n-array, expected return(s) (n)                 double
c        covb   : n-array, covariance asset vs.index (n)          double  
c        varb   : variance of index/benchmark                     double
c        rhob   : expected return of index/benchmark              double
c        neq    : number equality constraints                    integer
c        nin    : number inequality constraints                  integer
c        ccst   : matrix of linear constraints (n*(neq+nin))      double
c        bcst   : vector of linear constraints (neq+nin)          double
c        cinf   : lower bounds (n)                                double
c        csup   : upper bounds (n)                                double
c        temax  : tracking error constraint (>0)                  double
c        maxdic : maximum iteration (dichotomy)                  integer
c        epsdic : precision (dichotomy)                           double
c        epskat : Kato sensibility parameter                      double
c
c     WORKSPACE 
c        iwork  : vector ( 16*n + neq + 2*nin + 1 )              integer 
c        dwork  : matrix                                          double
c                 n*(7*n + neq + nin + 42) 
c                 + 2*neq + 3*nin + 4
c
c     OUTPUT 
c        wopt   : optimal portfolio vector (n)                   double
c        delopt : 
c        optte  : optimal Tracking Error                         double
c        info   : diagnostic argument                           integer
c
c     CALL   
c        WZERO  : put zero at the values of a vector less than eps
c                 and increase the greatest element with residual
c        YV     : copy a vector in a vector
c        XV     : computing scalar product of two vectors = V1'*V2
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, neq, nin, maxiter, info
      DOUBLE PRECISION rhob, varb, temax, epsilon, optte, delopt
      DOUBLE PRECISION cov(*), rho(*), cinf(*), csup(*), ccst(*), 
     &                 bcst(*), wopt(*), covb(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, niter, pdwopt, pite, pdte, piwq, pdwq, pccs, 
     &        pbcs, plin, plagr, neqtot, nintot, infoq, infote, infotmp
      DOUBLE PRECISION emret, tesav, delsav, epste, eps, delmin, delmax,
     &                 rdtmax 
      PARAMETER ( eps = 1.0E-30 )
c
c     external subroutines
      EXTERNAL WZERO, YV, XV, ctit, qp
c     
c     intrinsic functions
      INTRINSIC abs       
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwq = 1
c     piwq  : pointer for internal workspaces of QUAPRO
c             needs ( 3*n + neq + 2*nin + 1 )
      pite = piwq + ( 3*n + neq + 2*nin + 1 ) 
c     pite  : pointer for EXATER, needs*( 13*n )    
c
c     Total size of iwork array = ( 16*n + neq + 2*nin + 1 )
c     
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwopt = 1
c     pdwopt : pointer for wopt save, need ( n )
      pccs = pdwopt + ( n ) 
c     pccs   : pointer for matrix constraint in CTIT
c              needs ( n *( neq + nin + 1)) 
      pbcs = pccs + ( n *( neq + nin + 1)) 
c     pbcs   : pointer for vector constraint in CTIT
c              needs ( neq + nin + 1 )
      plin = pbcs + ( neq + nin + 1 )
c     plin   : pointer for linear part in CTIT, needs ( n )    
      pdwq   = plin + ( n )
c     pdwq   : pointer for workspace of QUAPRO
c              needs ( n*n + 6*n + nin ) 
      plagr  = pdwq + ( n*n + 6*n + nin )
c     plagr  : pointer for Lagrangian QUARPO        
c              needs ( n + neq + nin )
      pdte   = pdwq + ( n + neq + nin )
c     pdte   : pointer for EXATER needs ( n*(6*n+32) + 3 ) 
c
c     Total size of dwork = n*(7*n + neq + nin + 42 ) 
c                           + 2*neq + 3*nin + 4
c-----------------------------------------------------------------------
c
c     computes max of returns
      rdtmax = rho(1)
      DO i=2,n
         rdtmax = MAX( rdtmax, rho(i) )
      ENDDO
      delmax  = rdtmax - rhob
c
c     initializations
      infotmp = 0
      tesav   = 0
      delsav  = 0
      delopt  = -delmax
      delmin  = -delmax
      epste   = 2.*epsilon
      optte   = 0
      niter   = 0
c
c     stop test: (niter >= maxiter)&(|optTE - maxTE| <= epsilon)      
      DO WHILE ( ((niter+1).le. maxiter).and.(epste .gt. epsilon) )
c
c      increase number of iteration
       niter = niter + 1
c
c      constraints matrix and linear part
       CALL ctit ( n, rhob, rho, covb,
     &            delopt, neq, nin, ccst, bcst,
     &            neqtot, nintot, dwork(pccs), dwork(pbcs), dwork(plin))
c
c      quadratic solver (QP)
c      min [w*Q*w - p*w] s.t. C*w <= b, cinf <= w <= csup
       CALL qp ( n, cov, dwork(plin), dwork(pccs),
     &           dwork(pbcs), cinf, csup, neqtot, nintot,
     &           iwork(piwq), dwork(pdwq), dwork(plagr), wopt,
     &           infoq )
c
c       if wopt(i)< myzero -> wopt(i) = 0
c                           -> increase max[wopt(i)]
        CALL WZERO ( n, wopt, eps )
c
c       optimal ex-ante mean return = w'*rho
        CALL XV ( n, wopt, rho, emret )
c
c       optimal over performance
        delopt = emret - rhob
c
c       optimal ex-ante Tracking Error
        CALL EXATER ( n, wopt, cov, covb, varb, dwork(pdte),
     &                optte, infote)
c
c       stop test: epste = |optTE - maxTE|
        epste = abs(optte - temax)
c
c       error management
        IF (infoq .LT. 0) THEN         
           IF (infoq .ne. 1001) THEN
              RETURN
           ELSE
              infotmp = 1001
           ENDIF
        ELSE
           infotmp = 0
           delsav  = delopt
           tesav   = optte
           CALL YV ( n, wopt, dwork(pdwopt) )
        ENDIF
c
c       dichotomy (case niter=1 and niter>1)
        IF (niter .EQ. 1) THEN 
           IF (temax .ge. optte) THEN
               delmin = delopt - epsilon
           ELSE
               info = 107
               RETURN 
           ENDIF
        ENDIF
c   
        IF (niter .GT. 1) THEN
           IF (temax .GT. optte) THEN 
               delmin = delopt - epsilon
           ELSE
               delmax = delopt + epsilon           
           ENDIF
        ENDIF
c         
c       dichotomy -> 
        delopt = (delmax + delmin)/2.0     
c
c     end while
      ENDDO
c      
      IF (infotmp .eq. 0) THEN
           info   = infotmp
           delopt = delsav
           optte  = tesav
           CALL YV ( n, dwork(pdwopt), wopt )
      ENDIF
c      
      IF ((niter + 1) .gt. maxiter) info = 106 
c      
      RETURN
      END
