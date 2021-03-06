c=======================================================================
c
c     SOLSCOR
c
c     Computing the solution of SDLS optimization from the dual solution
c     specific to obtain a correlation matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE solscor ( nmat, cmat, ydual, iwork, dwork,
     &                     xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix (nmat*nmat) to optimize
c                     vectorized  vector(nmat*(nmat+1)/2)         double
c            mct    : number of constraints                      integer
c            ydual  : dual solution                   vector(mct) double
c
c     WORKSPACE 
c            iwork  : vector (12*nmat)                           integer
c            dwork  : vector (nmat*(nmat+1)/2 + nmat*(3*nmat+27)) double
c
c     OUTPUT 
c            xmat   : symmetric matrix solution
c                     vectorized  vector(nmat*(nmat+1)/2)         double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        CVMDS    : converting a vector in vectorized symmetric
c                   diagonal matrix
c        SMS      : computing the sum of 2 vectorized symmetric matrices
c        PROJECTS : projection of a symmetric matrix
c                   on the cone of semidefinite matrices
c                   version with vector (n*(n+1)/2) in input/output
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nmat, info
      DOUBLE PRECISION cmat(*), ydual(*), xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER ssm, piw, pmat, pdw
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      ssm = (nmat*(nmat+1)/2)
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw    = 1
c     piw    : pointer for workspace used by PROJECTS, so nmat more
c     pinext : pointer for the next iwork array
c              pinext = piw + 12*nmat
c              
c     Total size of iwork array (12*nmat)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat  : pointer for a work vectorized symmetric matrix
c             so ssm = (nmat*(nmat+1)/2) more
      pdw = pmat + ssm
c     pdw   : pointer for workspace used by PROJECTS
c             so (nmat*(3*nmat+27))  more
c     pnext :  pointer for the next work array
c              pnext = pdw + (nmat*(3*nmat+27)
c                    = 1 + ssm + (nmat*(3*nmat+27)) 
c
c     Total size of work array
c                    = (nmat*(nmat+1)/2) + (nmat*(3*nmat+27))
c
c-----------------------------------------------------------------------
c
c     computing the adjoint : symmetric matrix sum of Ai*y(i)
c     for correlation : gives a diagonal matrix with y 
      CALL CVMDS ( nmat, ydual, dwork(pmat) )
c
c     computing the sum C + adjoint = ( C+Aiy(i) )
      CALL SMS ( nmat, cmat, dwork(pmat), dwork(pmat) )
c
c     projection x = pK( C+Aiy(i) )
      CALL PROJECTS ( nmat, dwork(pmat), iwork(piw), dwork(pdw),
     &                xmat, info  )   
c
      RETURN
      END
