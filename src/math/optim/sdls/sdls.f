c=======================================================================
c
c     subroutine  SDLS                                       
c
c     Semi-Definite Least-Square optimization
c     with M >= alpha >= 0.0
c     ( the constraints matrices Ai may be non-symmetric )
c
c-----------------------------------------------------------------------
      SUBROUTINE sdls ( n, cmat, mcte, acte, bcte, epsbfg, alpha,
     &                  iwork, dwork, xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                     integer
c       cmat   : symmetric matrix (n*n) to optimize               double
c       mcte   : number of equal constraints                     integer
c       acte   : mcte Ai symmetric matrices (n*n) of
c                     equal constraints, vector(mcte*n*n)         double
c       bcte   : value of equal constraints, vector(mcte)         double
c       epsbfg : BFGS precision stop test                         double
c       alpha  : alpha  ( 0.0 =<alpha )                           double
c
c     WORKSPACE 
c       iwork  : 2*mcte + 12*n + 3                               integer
c       dwork  :    mcte*(mcte + 25)/2 
c                + (2*mcte + 5)*(n*(n + 1)/2) 
c                +  n*(4*n + 27)                                  double
c     OUTPUT 
c       xmat   : symmetric matrix (n*n)                           double
c       info   : diagnostic argument                             integer
c
c     CALL   
c        SEMDX : computing the sum of the diagonal elements of a
c                matrix with a scalar : B = A + x*I
c                ( A square matrix(n*n), x scalar,
c                    gives B square matrix(n*n) )
c        YM    : copy a vectorized matrix in a vectorized matrix
c        TM    : computing the trace of a full matrix
c        YV    : copy a vector in a vector
c        CMCMS : converting a vector of full square matrix
c                in vector of symmetric matrix
c        ZMCMS : symmetrises a vectorized full matrix
c                and gives a vectorized symmetric matrix A=(A+A')/2
c        SDLSS : semi-Definite Least Square optimization
c                ( the constraints matrices Ai are symmetric )
c        CMSMC : converting a vector of symmetric matrix
c                in vector of full square matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, mcte, info
      DOUBLE PRECISION epsbfg, alpha
      DOUBLE PRECISION cmat(*), acte(*), bcte(*), xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, kb, kai, ksai, sm, ssm, pdbv, pdcv, pdcas, pdxs, pdsai,
     &        piw, pdw, pdev, pdvs
      DOUBLE PRECISION EPS, trace
      PARAMETER ( EPS = 1.E-15 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      sm   = n*n
      ssm  = n*(n + 1)/2
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw   = 1
c     piw   : pointer for SDLSS workspaces ( 2*mcte + 12*n + 3 )
c
c     Total size of iwork array : ( 2*mcte + 12*n + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdbv  = 1
c     pdbv  : pointer for changed variable vector b(mcte)
      pdcv  = pdbv + (mcte)
c     pdcv  : pointer for changed variable matrix C(n*n)
      pdcas = pdcv + ( n*n )
c     pdcas : pointer for the vectorized symmetric matrix C,
c             so ssm = (n*(n+1)/2) more
      pdxs  = pdcas + ssm
c     pdxs  : pointer for the vectorized symmetric matrix X,
c             so ssm = (n*(n+1)/2) more
      pdsai = pdxs + ssm
c     pdsai : pointer for the mcte vectorized symmetric matrix
c             Ai(n*(n+1)/2), so mcte*(n*(n+1)/2) more

      pdev  = 1
c     pdev : pointer for eigenvalues needs n
c      so n 
      pdvs = pdev + n
c     pdvs : pointer for eignevector needs n*n
c     so n*n + n 
c     which is obviously lower than pdsai +ssm 
      
      pdw   = pdsai + mcte*ssm
c     pdw   : pointer for local workspaces,
c             SDLSS needs (  (mcte*(mcte+23)/2)
c                         +  (mcte+3)*(n*(n+1)/2) + (n*(3*n+27))  )
c
c             eigenvd needs n*(n+6)
c
c
c     Total size of dwork array : (mcte) + (n*n) + ssm + ssm + mcte*ssm
c                               + (mcte*(mcte+23)/2)
c                               + (mcte+3)*(n*(n+1)/2)
c                               + (n*(3*n+27)) )
c
c                         with ssm = n*(n+1)/2
c
c     = mcte*(mcte+25)/2 + (2*mcte+5)*(n*(n+1)/2) + n*(4*n+27) 
c
c-----------------------------------------------------------------------
c
c     changing variable CA = C - alpha*I if alpha > epsilon
      IF ( alpha .GT. EPS ) THEN
         CALL SEMDX ( n, cmat, -alpha, dwork(pdcv) )
      ELSE
         CALL YM ( n, n, cmat, dwork(pdcv) )
      ENDIF
c
c     changing variable bai = bi - alpha*Tr(Ai) if alpha > epsilon
      IF (mcte .NE. 0) THEN
         IF (alpha .GT. EPS) THEN
            kai = 1
            kb  = pdbv
            DO i = 1,mcte
               CALL TM ( n, acte(kai), trace )
               dwork(kb) = bcte(i) - alpha*trace
               kai       = kai + sm
               kb        = kb + 1
            ENDDO
         ELSE
            CALL YV ( mcte, bcte, dwork(pdbv) )
         ENDIF
      ENDIF
c
c     converting CA : full matrix in vector of symmetric matrix
      CALL CMCMS ( n, dwork(pdcv), dwork(pdcas) )
c
c     symmetrising Ai matrices and converting in vector of symmetric
c     matrices
      IF (mcte .NE. 0) THEN
         kai  = 1
         ksai = pdsai
         DO i = 1,mcte
            CALL ZMCMS ( n, acte(kai), dwork(ksai) )
            kai  = kai + sm
            ksai = ksai + ssm
         ENDDO
      ENDIF
c
c     SDLS optimization (cf. sdlss.f)
      CALL sdlss ( n, dwork(pdcas), mcte, dwork(pdsai), dwork(pdbv),
     &             epsbfg, iwork(piw), dwork(pdw),
     &             dwork(pdxs), info )
c
c     converting a vector of symmetric matrix in vector of full matrix
      CALL CMSMC ( n, dwork(pdxs), dwork(pdcv) )
c
c     chgt. variable X = Z + alpha*I if alpha > epsilon
      IF (alpha .GT. EPS) THEN
         CALL SEMDX ( n, dwork(pdcv), alpha, xmat )
      ELSE
         CALL YM ( n, n, dwork(pdcv), xmat )
      ENDIF
c
c     test equalities constraints
      IF ((info .GE. 0).AND.(mcte .GE. 1)) THEN
         CALL checksdlseq ( n, xmat, mcte, acte, bcte,
     &                      dwork(pdw), info)
      ENDIF
      RETURN
      END
