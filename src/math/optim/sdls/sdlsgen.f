c=======================================================================
c
c     subroutine  SDLSGEN                                    
c
c     Semi-Definite Least Square optimization
c     general version with equal and inequal constraints
c     with M >= alpha >= 0.0
c     ( the constraints matrices Ai and Bi may be non-symmetric )
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlsgen ( n, cmat, mcte, acte, bcten,
     &                     mcti, acti, lcti, ucti, epsbfg, alpha,
     &                     iwork, dwork, xmat, info )
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
c        SDLS  : semi-Definite Least Square optimization
c                  general version with equal and inequal constraints
c                  ( the constraints matrices Ai are symmetric )
c        SDLSG  : semi-Definite Least Square optimization
c                  general version with equal and inequal constraints
c                  ( the constraints matrices Ai are symmetric )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, mcte, mcti, info
      DOUBLE PRECISION epsbfg, alpha
      DOUBLE PRECISION cmat(*), xmat(*), acte(*), bcten(*), acti(*),
     &                 lcti(*), ucti(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pdw
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c
c     Total size of iwork array : (2*mcte+4*mcti+12*n+5)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdw  = 1
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
c-----------------------------------------------------------------------
c
      IF (mcti .EQ. 0) THEN
        CALL sdls ( n, cmat, mcte, acte, bcten, epsbfg, alpha,
     &              iwork(piw), dwork(pdw), xmat, info )
      ELSE
        CALL sdlsg ( n, cmat, mcte, acte, bcten,
     &               mcti, acti, lcti, ucti, epsbfg, alpha,
     &               iwork(piw), dwork(pdw), xmat, info )
      ENDIF
c
      RETURN
      END
