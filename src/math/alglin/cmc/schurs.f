c=======================================================================
c
c     subroutine SCHURS                                      
c
c     Matrix Schur Decomposition of a symmetric matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE schurs ( n, X, iwork, dwork, eigval, vpmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : input matrix (n*n)                          double
c
c     WORKSPACE 
c            iwork  : vector ( 12*n )                            integer
c            dwork  : vector ( n*(n+26) )                         double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            vpmat  : eigenvectors (n*n)                          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            YM     : copy a vectorized matrix in a vectorized matrix
c            DSYEVR : Schur factorization for symmetric matrices (cf. LAPACK)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION eigval(*), X(*), vpmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pisup, pmat, pdw, ldwork, liwork, ild, ili, il, iu, 
     &        msy
      DOUBLE PRECISION vl, vu, abstol
      PARAMETER ( abstol = 0.0, ili = 10, ild = 26 )
      PARAMETER ( il = 1, vl = 0.0, vu = 0.0 )
c
c     external subroutines
      EXTERNAL YM, DSYEVR
c
c-----------------------------------------------------------------------
c
c     initialisations
      info = 0 
      liwork = ili*n
      ldwork = ild*n
      iu = n
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for work matrix of DSYEVR who needs  (n*ili)
c              so n*ili = n*10 more
      pisup = piw + ( 10*n )
c     pisup : pointer for work matrix of DSYEVR who needs  (2*n) so 2*n more
c     pinext = pisup + ( 2*n )
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( 12*n )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat   : pointer for work matrix, so (n*n) more
      pdw = pmat + ( n*n )
c     pdw    : pointer for work matrix of DSYEVR ( n*ild ),
c              so n*ild = n*26 more
c     pdnext = pdw + ( 26*n ) 
c     pdnext  : pointer of the next dwork array
c
c     Total size of dwork array = ( n*(n + 26) )
c
c-----------------------------------------------------------------------
c
c     saves initial matrix ( altered by DSYEVR )
      CALL YM ( n, n, X, dwork(pmat) )
c
c     DSYEVR (cf. LAPACK) 
c     DSYEVR computes selected eigenvalues and, optionally, eigenvectors
c     of a real symmetric matrix T.  Eigenvalues and eigenvectors can be
c     selected by specifying either a range of values or a range of
c     indices for the desired eigenvalues
      CALL DSYEVR ( 'V', 'A', 'U', n, dwork(pmat), n,
     &              vl, vu, il, iu, abstol, msy, eigval,
     &              vpmat, n, iwork(pisup),
     &              dwork(pdw), ldwork, iwork(piw), liwork, info )
c
      RETURN
      END

