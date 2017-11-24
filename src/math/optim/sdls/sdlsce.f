c=======================================================================
c
c     subroutine  SDLSCE
c
c     Semi-Definite Least Square optimization
c     with equal constraints and M >= alpha >= 0.0
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlsce ( nmat, cmat, imatce, matce, epsbfg, alpha,
     &                    iwork, dwork, xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            imatce : symmetric matrix(nmat*nmat)                integer
c                     equal constraints indicators
c                     if imatce(i,j) = 1 : xmat(i,j) = cmat(i,j)
c                     if imatce(i,j) = 2 : xmat(i,j) = matce(i,j)
c                     only the lower part is used
c            matce  : symmetric matrix(nmat*nmat)                 double
c                     equal constraint value if imatce(i,j) = 2
c                     only the lower part is used
c            epsbfg : BFGS precision stop test                    double
c            alpha  : alpha  ( 0.0 =<alpha )                      double
c
c     WORKSPACE 
c            iwork  : workspace                                  integer
c                     vector( 2*mcte + 12*nmat + 3 )
c            dwork  : vector (   mcte*(nmat*nmat)
c                              + mcte*(mcte+27)/2
c                              + (2*mcte+5)*(nmat*(nmat+1)/2)
c                              + (nmat*(4*nmat+27))  )            double
c                     mcte is the number of imatce lower part
c                     elements different of zero
c     OUTPUT 
c            xmat   : symmetric matrix solution (nmat*nmat)       double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        CSDLSCE : computing the equal constraints matrices for SDLSCE
c        SDLS    : semi-Definite Least Square optimization
c                  with M >= alpha >= 0.0
c                  ( the constraints matrices Ai may be non-symmetric )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, info, imatce(nmat,*)
      DOUBLE PRECISION epsbfg, alpha, cmat(*), matce(*), xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mcte, i, j, pdb, pda, piw, pdw
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     number of equalities constraints
      mcte = 0
      DO i = 1,nmat
         DO j = 1,i
            IF ( imatce(i,j).ne.0 ) mcte = mcte + 1
         ENDDO
      ENDDO
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw   = 1
c     piw   : pointer for local workspaces,
c             SDLS needs ( 2*mcte + 12*nmat + 3 )
c
c     Total size of iwork array : ( 2*mcte + 12*nmat + 3 )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdb   = 1
c     pdb   : pointer for vector b(mcte)
      pda   = pdb + (mcte)
c     pda   : pointer for the mcte matrix(nmat*nmat),
c             so mcte*(nmat*nmat) more
      pdw   = pda + mcte*(nmat*nmat)
c     pdw   : pointer for local workspaces,
c             SDLS needs (  (mcte*(mcte+25)/2)
c                         + (2*mcte+5)*(nmat*(nmat+1)/2)
c                         + (nmat*(4*nmat+27))
c
c     Total size of dwork array : (mcte) + mcte*(nmat*nmat)
c                               + (mcte*(mcte+25)/2)
c                               + (2*mcte+5)*(nmat*(nmat+1)/2)
c                               + (nmat*(4*nmat+27))
c
c     =  mcte*(nmat*nmat)
c      + mcte*(mcte+27)/2
c      + (2*mcte+5)*(nmat*(nmat+1)/2)
c      + (nmat*(4*nmat+27))
c
c-----------------------------------------------------------------------
c
c     computing Ai constraints matrices (cf. utsdls.f)
      CALL csdlsce ( nmat, cmat, imatce, matce,
     &               mcte, dwork(pda), dwork(pdb) )
c
c     SDLS optimization
      CALL sdls ( nmat, cmat, mcte, dwork(pda), dwork(pdb),
     &            epsbfg, alpha, iwork(piw), dwork(pdw),
     &            xmat, info )
c
      RETURN
      END
