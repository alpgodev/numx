c=======================================================================
c
c     subroutine SCHUR                                       
c
c     Matrix Schur Decomposition
c
c     If A is a complex square matrix of dimention n, then there exists 
c     a unitary matrix Q such that
c
c     Q'AQ = T = D + N
c
c     where ' is the conjugate transpose, 
c     T is an upper triangular matrix,
c     D = Diag(l1,...,lp) (the li are eigenvalues of A), 
c     and N is strictly upper triangular matrix. 
c
c     Furthermore, Q can be chosen such that the eigenvalues li appear 
c     in any order along the diagonal.
c
c-----------------------------------------------------------------------
      SUBROUTINE schur ( n, X, iwork, dwork, eigval, vpmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : input matrix (n*n)                          double
c
c     WORKSPACE 
c            iwork  : vector (n)                                 integer
c            dwork  : vector ( n*(n+6) )                          double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            vpmat  : eigenvectors (n*n)                          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            YM     : copy a vectorized matrix in a vectorized matrix
c            DGEES  : Schur factorization (cf. LAPACK)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION eigval(*), X(*), vpmat(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      LOGICAL select
      INTEGER sdim, piw, pmat, pdw, pip, lwork, ilwork
      PARAMETER ( ilwork = 5 )
c
c     external subroutines
      EXTERNAL YM, DGEES
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     lwork  : size of workspace for DGEES 
c              ( min 3*N , so ilwork = min 3 , more is best )
      lwork = n*ilwork
      sdim = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for DGEES  who needs (n)
c     pinext = piw + n
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( n )
c
c     pointers for double precision work space  : work(nasmax,*)
c     ----------------------------------------------------------
      pmat = 1
c     pmat   : pointer for work  matrix so (n*n) more
      pdw = pmat + ( n*n ) 
c     pdw    : pointer for work matrix of DGEES so lwork=(n*5) more
      pip = pdw + lwork
c     pip    : pointer for imaginary part vector so (n) more
c     pdnext = pip + n
c     pdnext  : pointer of the next dwork array
c
c     Total size of work array = ( n*(n + 6) )
c
c-----------------------------------------------------------------------
c
c     saves matrix initial ( altered by DGEES )
      CALL YM ( n, n, X, dwork(pmat) )
c
c     DGEES (cf. LAPACK)
c     DGEES computes for an N-by-N real  nonsymmetric  matrix  A,
c     the eigenvalues, the real Schur form T, and, optionally, the
c     matrix of Schur vectors Z
      CALL DGEES( 'V', 'N', select, n, dwork(pmat), n,
     &            sdim, eigval, dwork(pip), vpmat,
     &            n, dwork(pdw), lwork, iwork(piw), info )
c
      RETURN
      END

