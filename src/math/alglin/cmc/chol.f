c=======================================================================
c
c     subroutine CHOL
c
c     Choleski factorization of a symmetric and definite matrix:
c
c     X = L*L', where L is an lower triangular matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE chol ( n, X, dwork, L, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : size of matrix X                           integer
c            X      : symmetric matrix (n*n)                      double
c
c     WORKSPACE 
c            dwork  : vector (n*n)                                double
c
c     OUTPUT 
c            L      : lower triangular matrix (n*n)               double                   
c            info   : diagnostic argument                        integer
c
c     CALL   
c            YM     : copy a vectorized matrix in a vectorized matrix
c            DPOTRF : Choleski factorization (cf. LAPACK)
c            CMCML  : converting a full square matrix (n*n)
c                    in low triangular matrix (n*n)
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
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pmat
c
c     external subroutines
      EXTERNAL YM, DPOTRF, CMCML     
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat  : pointer for vectorized symmetric matrix, so (n*n) more
c
c     Total size of dwork array = (n*n)
c
c-----------------------------------------------------------------------
c
c     saves matrix initial ( altered by DPOTRF )
      CALL YM ( n, n, x, dwork(pmat) )
c
c     Choleski factorization (cf. LAPACK)
      CALL DPOTRF ( 'L', n, dwork(pmat), n, info )
      IF (info .ne. 0) RETURN
c
c     conversion in low triangular matrix
      CALL CMCML ( n, dwork(pmat), L)
c
      RETURN
      END

