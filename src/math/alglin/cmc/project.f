c=======================================================================
c
c     subroutine PROJECT                                     
c
c     Symmetric matrix projection on the semidefinite matrix cone
c
c-----------------------------------------------------------------------
      SUBROUTINE project ( n, X, iwork, dwork, promat, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : matrix size                                integer
c            X      : symmetric matrix (n*n)                      double
c
c     WORKSPACE :
c            iwork  : vector ( 12*n )                            integer
c            dwork  : vector ( n*( 3*n + 27 ) )                   double
c
c     OUTPUT :
c            promat : projected matrix (n*n)                      double
c            info   : diagnostic argument                        integer
c
c     CALL   :
c        SCHURS  : SCHUR decomposition of a symmetric matrix
c        OMCDMCT : computing M*D*M'
c                  (M square matrix n*n, D vector n of diagonal matrix,
c                  gives square matrix n*n )
c        ZMCMC   : symmetrises a vectorized full matrix
c                  and gives a vectorized full square matrix A=(A+A')/2
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments 
      INTEGER n, info
      DOUBLE PRECISION X(*), promat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pmat, pev, pvs, pdw, pdm, i
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL schurs, OMCDMCT, ZMCMC
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for work space of SCHURS, so ( 12*n ) more
c
c     Total size of iwork array ( 12*n )
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pmat = 1
c     pmat   : pointer for work  matrix ( n*n ), so n*n more
      pev = pmat + ( n*n )
c     pev    : pointer eigenvalues vector ( n ), so n more
      pvs = pev + ( n )
c     pvs    : pointer eigenvectors matrix ( n*n ), so n*n more
      pdw = pvs + ( n*n )
c     pdw    : pointer for work space of SCHURS, so n*(n + 26) more
      pdm = pdw
c     pdm    : pointer for work matrix of OMCDMCT ( n*n ),
c              uses the same space than SCHURS, so nothing more
c
c              dwtot = 2*n*n + n + n*(n+26)
c                    = n*( 3*n + 27 )
c
c     Total size of work array = ( n*( 3*n + 27 ) )
c
c----------------------------------------------------------------------
c
c     Schur factorization
      CALL schurs ( n, X, iwork(piw), dwork(pdw),
     &              dwork(pev), dwork(pvs), info )
      IF (info .lt. 0) RETURN
c
c     projection on the cone of semidefinite matrices
      DO i = 1,n
         IF ( dwork(pev+i-1) .LT. EPS ) dwork(pev+i-1) = EPS
      ENDDO
c
c     reconstruction of the matrix = Z*T*Z'
      CALL OMCDMCT ( n, dwork(pvs), dwork(pev), dwork(pdm), dwork(pmat))
c
c     assures the symmetry of the matrix (computer precision)
      CALL ZMCMC ( n, dwork(pmat), promat )
c
      RETURN
      END

