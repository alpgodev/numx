c=======================================================================
c
c     subroutine  RCHO                                       
c
c     Robust Choleski factorization of a symmetric matrix (not necessary 
c     definite positive): 
c
c     X = L*L', where L is an lower triangular matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE rcho ( n, X, iwork, dwork, L, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : symmetric matrix (n*n)                      double
c
c     WORKSPACE 
c            iwork  : vector(12*n)                               integer
c            dwork  : vector(n*(4*n + 27))                        double
c
c     OUTPUT 
c            L      : lower triangular matrix (n*n)               double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c            PROJECTA : projection of a symmetric matrix on matrix cone >= alpha
c            CHOL     : Choleski factorization X = L*L' of a symmetric matrix 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION X(*), L(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pmat, pdw
      DOUBLE PRECISION alpha
      PARAMETER ( alpha = 1.E-30 )
c
c     external subroutines
      EXTERNAL projecta, chol
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for PROJECTA workspace who needs ( 12*n )
c
c     Total size of iwork array : ( 12*n )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat  : pointer for work matrix(n*n), so (n*n) more
      pdw = pmat + ( n*n )
c     pdw   : pointer for PROJECTA workspace who needs n*(3*n+27)
c                        and CHOL workspace who needs (n*n),
c             so n*(3*n+27) more
c
c     Total size of dwork array : (n*n) + n*(3*n+27)
c                                 = n*(4*n + 27)
c
c-----------------------------------------------------------------------
c
c     projection on the cone of definite matrices
      CALL projecta ( n, X, alpha, iwork(piw), dwork(pdw),
     &                dwork(pmat), info )
      IF (info .LT. 0) RETURN
c
c     Choleski factorization
      CALL chol ( n, dwork(pmat), dwork(pdw), L, info )
c
      RETURN
      END

