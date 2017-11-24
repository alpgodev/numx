c=======================================================================
c
c     subroutine  SDLSG
c
c     Semi-Definite Least Square optimization
c     general version with equal and inequal constraints
c     with M >= alpha >= 0.0
c     ( the constraints matrices Ai and Bi may be non-symmetric )
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlsg ( n, cmat, mcte, acte, bcte,
     &                   mcti, acti, lcti, ucti, epsbfg, alpha,
     &                   iwork, dwork, xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                     integer
c       cmat   : matrix (n*n) to optimize                         double
c       mcte   : number of equal constraints                     integer
c       acte   : mcte Ai symmetric matrices (n*n) of
c                equal constraints  vector(mcte*n*n)              double
c       bcte   : vector(mcte) of equal constraints                double
c       mcti   : number of inequal constraints                   integer
c       acti   : mcti Bi symmetric matrices(n*n) of
c                inequal constraints  vector(mcti*n*n)            double
c       lcti   : lower bounds of inequal constraints
c                vector(mcti)                                     double
c       ucti   : upper bounds of inequal constraints
c                vector(mcti)                                     double
c       epsbfg : BFGS precision stop test                         double
c       alpha  : alpha  ( 0.0 =<alpha )                           double
c
c     WORKSPACE 
c       iwork  : 2*mcte + 4*mcti + 12*n + 5                      integer
c       dwork  : vector                                           double
c                (   (mcte+2*mcti)*(mcte+2*mcti+25)/2
c                  + (2*mcte+2*mcti+5)*(n*(n+1)/2)
c                  + (n*(4*n+27))  )
c
c     OUTPUT 
c       xmat   : matrix solution vectorized vector(n*n)          double
c       info   : diagnostic argument                            integer
c
c     CALL   
c        SEMDX   : computing the sum of the diagonal elements of a
c                  matrix with a scalar : B = A + x*I
c                  ( A square matrix(n*n), x scalar,
c                    gives B square matrix(n*n) )
c        YM      : copy a vectorized matrix in a vectorized matrix
c        TM      : computing the trace of a full matrix
c        YV      : copy a vector in a vector
c        CMCMS   : converting a vector of full square matrix
c                  in vector of symmetric matrix
c        ZMCMS   : symmetrises a vectorized full matrix
c                  and gives a vectorized symmetric matrix A=(A+A')/2
c        SDLSSG  : semi-Definite Least Square optimization
c                  general version with equal and inequal constraints
c                  ( the constraints matrices Ai are symmetric )
c        CMSMC   : converting a vector of symmetric matrix
c                  in vector of full square matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simsdls
c
c     i/o arguments
      INTEGER n, mcte, mcti, info
      DOUBLE PRECISION epsbfg, alpha
      DOUBLE PRECISION cmat(*), xmat(*), acte(*), bcte(*), acti(*),
     &                 lcti(*), ucti(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, kb, kl, ku, ki, ksi, sm, ssm, pdbv, pdlv, pduv, pdcv,
     &        pdcas, pdxs, pdsai, pdsbi, piw, pdw
      DOUBLE PRECISION eps, trace
      PARAMETER ( eps = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      sm   = n*n
      ssm  = n*(n+1)/2
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for local workspaces,
c             SDLSS needs (2*mcte+4*mcti+12*n+5)
c
c     Total size of iwork array : (2*mcte+4*mcti+12*n+5)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdbv  = 1
c     pdbv  : pointer for changed variable vector b(mcte)
      pdlv  = pdbv + (mcte)
c     pdlv  : pointer for changed variable vector low(mcti)
      pduv  = pdlv + (mcti)
c     pduv  : pointer for changed variable vector up(mcti)
      pdcv  = pduv + (mcti)
c     pdcv  : pointer for changed variable matrix C(n*n)
      pdcas = pdcv + (n*n)
c     pdcas : pointer for the vectorized symmetric matrix C,
c             so ssm = (n*(n+1)/2) more
      pdxs  = pdcas + ssm
c     pdxs  : pointer for the vectorized symmetric matrix X,
c             so ssm = (n*(n+1)/2) more
      pdsai = pdxs + ssm
c     pdsai : pointer for the mcte vectorized symmetric matrix
c             Ai(n*(n+1)/2), so mcte*(n*(n+1)/2) more
      pdsbi = pdsai + mcte*ssm
c     pdsbi : pointer for the mcti vectorized symmetric matrix
c             Bi(n*(n+1)/2), so mcti*(n*(n+1)/2) more
      pdw   = pdsbi + mcti*ssm
c     pdw   : pointer for local workspaces,
c             SDLSSG needs (  (mcte+2*mcti)*((mcte+2*mcti+23)/2)
c                           + (mcte+mcti+3)*(n*(n+1)/2)
c                           + (n*(3*n+27))  )
c
c     Total size of dwork array : (mcte) + (mcti) + (mcti) + (nmat*nmat)
c                               + ssm + ssm + mcte*ssm + mcti*ssm
c                               + (mcte+2*mcti)*(mcte+2*mcti+23)/2
c                               + (mcte+mcti+3)*(n*(n+1)/2)
c                               + (n*(3*n+27))  )
c
c                         with ssm = n*(n+1)/2
c
c                         = (mcte+2*mcti)*(mcte+2*mcti+25)/2
c                         + (2*mcte+2*mcti+5)*(n*(n+1)/2)
c                         + (nmat*(4*n+27)) )
c
c-----------------------------------------------------------------------
c
c     chgt. variable CA = C - alpha*I (alpha > 0)
      IF (alpha .GT. eps) THEN
         CALL SEMDX ( n, cmat, -alpha, dwork(pdcv) )
      ELSE
         CALL YM ( n, n, cmat, dwork(pdcv) )
      ENDIF
c
c     chgt. variable bai = bi - alpha*Tr(Ai) (alpha > 0)
      IF (mcte .NE. 0) THEN
         IF (alpha .GT. eps) THEN
            ki = 1
            kb = pdbv
            DO i = 1,mcte
               CALL TM ( n, acte(ki), trace )
               dwork(kb) = bcte(i) - alpha*trace
               ki = ki + sm
               kb = kb + 1
            ENDDO
         ELSE
            CALL YV ( mcte, bcte, dwork(pdbv) )
         ENDIF
      ENDIF
c
c     chgt. variable lbi = li - alpha*Tr(Bi) (alpha > 0)
c     chgt. variable ubi = ui - alpha*Tr(Bi) (alpha > 0)
      IF (mcti .NE. 0) THEN
         IF (alpha .GT. eps) THEN
            ki = 1
            kl = pdlv
            ku = pduv
            DO i = 1,mcti
               CALL TM ( n, acti(ki), trace )
               dwork(kl) = lcti(i) - alpha*trace
               dwork(ku) = ucti(i) - alpha*trace
               ki = ki + sm
               kl = kl + 1
               ku = ku + 1
            ENDDO
         ELSE
            CALL YV ( mcti, lcti, dwork(pdlv) )
            CALL YV ( mcti, ucti, dwork(pduv) )
         ENDIF
      ENDIF
c
c     converting C : vector of full matrix in vector of symmetric matrix
      CALL CMCMS ( n, dwork(pdcv), dwork(pdcas) )
c
c     symmetrising Ai matrices and converting in vector of symmetric
c     matrices
      IF (mcte .NE. 0) THEN
         ki  = 1
         ksi = pdsai
         DO i = 1,mcte
            CALL ZMCMS ( n, acte(ki), dwork(ksi) )
            ki  = ki + sm
            ksi = ksi + ssm
         ENDDO
      ENDIF
c
c     symmetrising Bi matrices
      IF (mcti .NE. 0) THEN
         ki  = 1
         ksi = pdsbi
         DO i = 1,mcti
            CALL ZMCMS ( n, acti(ki), dwork(ksi) )
            ki  = ki + sm
            ksi = ksi + ssm
         ENDDO
      ENDIF
c
c     SDLS optimization
      CALL sdlssg ( n, dwork(pdcas), mcte, dwork(pdsai), dwork(pdbv),
     &              mcti, dwork(pdsbi), dwork(pdlv), dwork(pduv),
     &              epsbfg, iwork(piw), dwork(pdw),
     &              dwork(pdxs), info )
c
c     converting a vector of symmetric matrix in vector of full matrix
c     for solution X
      CALL CMSMC ( n, dwork(pdxs), dwork(pdcv) )
      
c
c     chgt. variable X = Z + alpha*I if alpha > epsilon
      IF (alpha .GT. eps) THEN
         CALL SEMDX ( n, dwork(pdcv), alpha, xmat )
      ELSE
         CALL YM ( n, n, dwork(pdcv), xmat )
      ENDIF
      IF (info .GE. 0) THEN
         CALL checksdls ( n, xmat, mcte, acte, bcte,
     &                    mcti, acti, lcti, ucti, dwork(pdw), info)
      ENDIF
c
      RETURN
      END
