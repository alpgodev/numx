c=======================================================================
c
c     subroutine SCHURSO                                     
c
c     Matrix Schur Decomposition of a symmetric matrix
c     and sorts eigenvalues (decrease order)
c
c-----------------------------------------------------------------------
      SUBROUTINE schurso ( n, X, iwork, dwork,
     &                     eigval, soreig, vpmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : input matrix (n*n)                          double
c
c     WORKSPACE 
c            iwork  : vector ( 12*n )                            integer
c            dwork  : vector ( n*(2*n+26) )                       double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            soreig : sorted eigenvalues (n)                      double
c            vpmat  : ordered eigenvectors (n*n)                  double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            SCHURS : Schur decomposition of a symmetric matrix
c            AVODI  : sorts a vector in decrease order with index in old vector
c            ACMI   : sorts the columns of a vectorized matrix as index
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION eigval(*), soreig(*), X(*), vpmat(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pvs, pdw
c
c     external subroutines
      EXTERNAL schurs, AVODI, ACMI
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw    : pointer for workspaces SCHURS needs  (12*n)
c              AVODI and ACMI need  (n), so (12*n) more
c     pinext = piw + 12*n
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( 12*n )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pvs = 1
c     pvs    : pointer eigenvectors matrix so (n*n) more
      pdw = pvs + ( n*n )
c     pdw    : pointer SCHURS who needs ( n*(n+26) ), so ( n*(n + 26) ) more
c     pdnext = pdw + ( n*(n+26) ) 
c     pdnext  : pointer of the next dwork array
c
c     Total size of dwork array = ( n*(2*n + 26) )
c
c-----------------------------------------------------------------------
c
c     Schur factorization
      CALL schurs ( n, X, iwork(piw), dwork(pdw),
     &              eigval, dwork(pvs), info )
      IF (info .lt. 0) RETURN
c
c     sorts the eigenvalues vector (decrease index order)
c     iwork is index of vector elements in the old vector
      CALL AVODI ( n, eigval, soreig, iwork(piw) )
c
c     sorts the columns of a vs (eigenvectors) as index
c     iwork is index of vector elements in the old vector
      CALL ACMI ( n, n, dwork(pvs), iwork(piw), vpmat )
c
      RETURN
      END

