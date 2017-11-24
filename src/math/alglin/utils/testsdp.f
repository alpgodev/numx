c=======================================================================
c
c     subroutine TESTSDP                                     
c
c     Matrix test if sdp (0) or not definite positive (1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE TESTSDP ( n, X, iwork, dwork, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : input matrix (n*n)                          double
c
c     WORKSPACE 
c            iwork  : vector (n)                                 integer
c            dwork  : vector ( n*(2*n+7) )                          double
c
c     OUTPUT 
c            info   : diagnostic argument                        integer
c                     =0, X is sdp 
c                     =1, X is not sdp
c
c     CALL   
c            DGEES  : Schur factorization (cf. LAPACK)
c            EVMIN  : minimum value of a vector (utils/utmat.f)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION X(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwk, pdwk, pdwev, pdwmat
      DOUBLE PRECISION minvp, EPS
      PARAMETER (EPS = 1.E-12)
c
c     external subroutines
      EXTERNAL schur, EVMIN
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info  = 0
      minvp = 0.0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwk = 1
c     piw   : pointer for SCHUR  who needs (n)
c
c     Total size of iwork array ( n )
c
c     pointers for double precision work space  : dwork
c     ----------------------------------------------------------
      pdwk = 1
c     pdwk   : pointer for SCHUR  so (n*(n + 6) 
      pdwev = pdwk + ( n*( n + 6 ) ) 
c     pdwev  : pointer for eigenvalues (n)
      pdwmat = pdwev + n
c     pdwmat : pointer for matrix (n*n) 
c
c     Total size of work array = ( n*(2*n + 7) )
c
c-----------------------------------------------------------------------
c
c     computation of the eigenvalues 
      CALL schur ( n, X, iwork(piwk), dwork(pdwk),
     &             dwork(pdwev), dwork(pdwmat), info )
c
c     minimum eigenvalue
      CALL EVMIN ( n, dwork(pdwev), minvp )
c
c     if spec(X) < 0 --> exit
      IF (minvp .LT. EPS) THEN
         info = 1
         RETURN
      ENDIF     
c
      RETURN
      END

