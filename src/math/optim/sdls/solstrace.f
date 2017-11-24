c=======================================================================
c
c     subroutine SOLSTRACE
c
c     Computing the solution of SDLS optimization from the dual solution
c     specific to Trace constraint
c
c-----------------------------------------------------------------------
      SUBROUTINE solstrace ( n, C, ydual, iwork, dwork, X, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                        integer
c       C      : matrix (n*n)                                        double
c       ydual  : dual solution                                       double
c
c     WORKSPACE 
c       iwork  : 12*n                                               integer
c       dwork  : n*(4*n+27)                                          double
c
c     OUTPUT 
c       X      : matrix solution (n*n)                               double
c       info   : diagnostic argument                                integer
c
c     CALL   
c       SM      : computing the sum of 2 matrices
c       PROJECT : projection of a symmetric matrix on SDP cone                 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, info
      DOUBLE PRECISION C(*), X(*), ydual
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pmat, pdw
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw    = 1
c     piw    : pointer for workspace used by PROJECT (n) 
c              
c     Total size of iwork array : 12*n
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat  : pointer for a work matrix (n*n)
      pdw = pmat + (n*n)
c     pdw   : pointer for workspace used by PROJECT (n*(3*n+27))
c
c     Total size of work array : n*(4*n+27)
c
c-----------------------------------------------------------------------
c
c     computing the sum: C + y*Id
      CALL SEMDX ( n, C, ydual, dwork(pmat) )
c      
c     projection X = pK( C + y*Id )
      CALL PROJECT ( n, dwork(pmat), iwork(piw), dwork(pdw), X, info )
c
      RETURN
      END
