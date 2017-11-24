c=======================================================================
c
c     subroutine SCHURO                                      
c
c     Matrix Schur Decomposition and sorts eigenvalues (decrease order)
c     cf. SCHUR subroutine 
c
c-----------------------------------------------------------------------
      SUBROUTINE schuro ( n, X, iwork, dwork,
     &                    eigval, soreig, vpmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : matrix (n*n)                                double
c
c     WORKSPACE 
c            iwork  : vector (n)                                 integer
c            dwork  : vector ( n*(2*n+6) )                        double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            soreig : sorted eigenvalues (n)                      double
c            vpmat  : eigenvectors (n*n)                          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            SCHUR  : matrix Schur decomposition
c            AVODI  : sorts a vector (decrease order) with index in old vector
c            ACMI   : sorts the columns of a vectorized matrix as index
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION eigval(*), soreig(*), X(*), vpmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pvs, pdw
c
c     external subroutines
      EXTERNAL schur, AVODI, ACMI
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for SCHUR, AVODI and ACMI  who need (n)
c     pinext = piw + ( n )
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( n )
c
c     pointers for double precision work space  : work(nasmax,*)
c     ----------------------------------------------------------
      pvs = 1
c     pvs    : pointer eigenvectors matrix so (n*n) more
      pdw = pvs + ( n*n )
c     pdw    : pointer SCHUR who needs ( n*(n+6) ), so ( n*(n+6) ) more
c     pdnext  : pointer of the next dwork array
c
c     Total size of work array = ( n*(2*n + 6) )
c
c-----------------------------------------------------------------------
c
c     Schur factorization
      CALL schur ( n, X, iwork(piw), dwork(pdw),
     &             eigval, dwork(pvs), info )
      IF (info .LT. 0) RETURN
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

