c=======================================================================
c
c     subroutine SCHURPCA
c
c     Matrix Schur Decomposition and sorts eigenvalues (decrease order)
c     cf. SCHUR subroutine 
c
c-----------------------------------------------------------------------
      SUBROUTINE schurpca ( n, M, iwork, dwork, eigval, eigvect, info)
c     
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            M      : matrix (n*n)                                double
c
c     WORKSPACE 
c            iwork  : vector (n)                                 integer
c            dwork  : vector ( n*(2*n+6) )                        double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            eigvect  : eigenvectors (n*n)                          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            SCHURO  : matrix Schur decomposition
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION M(*), eigval(*), eigvect(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pdunsorted, pdw
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for SCHURO who needs (n)
c
c     Total size of iwork array ( n )
c
c     pointers for double precision work space  : work(nasmax,*)
c     ----------------------------------------------------------
      pdunsorted = 1
c     pvs    : pointer unsorted vector so n more
      pdw = pdunsorted +  n
c     pdw    : pointer SCHUR0 who needs n*(2*n+6), so n*(2*n+6) more
c     pdnext  : pointer of the next dwork array
c
c     Total size of work array = n*(2*n + 7)
c
c-----------------------------------------------------------------------
c
c     Schur factorization
      CALL schuro ( n, M, iwork(piw), dwork(pdw), dwork(pdunsorted),
     &              eigval, eigvect, info )
c
      RETURN
      END
