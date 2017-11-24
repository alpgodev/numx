c=======================================================================
c
c     subroutine KATVAR
c
c     Kato decomposition of a symmetric matrix, maximum Kato sensibility 
c     and number factors for a given variance level
c
c-----------------------------------------------------------------------
      SUBROUTINE katvar ( n, X, vargiv, epsin, indpca,
     &                    iwork, dwork,
     &                    varopt, nbfact, epsout, soreig, vpmat,
     &                    nbkgrp, grpkat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : input matrix (n*n)                          double
c            vargiv : given variance (0.0 =< vargiv = <1.0)       double
c            epsin  : input Kato-sensitivity                      double
c            indpca : indicator (21, 22, 23)                     integer
c
c     WORKSPACE 
c            iwork  : vector (12*n)                              integer
c            dwork  : matrix ( n*(2*n + 27) )                     double
c
c     OUTPUT 
c            varopt : optimal variance                            double
c            nbfact : number of factors                          integer
c            epsout : maximum Kato sensibility                    double
c            soreig : sorted eigenvalues (n)                      double
c            vpmat  : ordered eigenvectors (n*n)                  double
c            nbkgrp : number of Kato groups                      integer
c            grpkat : Kato groups (n) output size (nbkgrp)       integer
c            info   : diagnostic argument                        integer
c
c     CALL   
c           SCHURSO : SCHUR decomposition of a symmetric matrix
c                     sorts eigenvalues and eigenvectors (decrease order)
c           GREPS   : computes groups for a sensibility epsilon
c
c     METHOD 
c            for Kato groups vector    
c            for an example of eigenvalues vector ( size = 7 ) :   
c            soreig = [ 5.2, 5.1, 3.5, 3.4, 3.3, 1.7 1.3 ]    
c            and eps = 0.15, the routine computes    
c            nbkgrp = 4    
c            grpvec = [ 2, 3, 1, 1 ]    
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, nbfact, nbkgrp, indpca, info
      INTEGER grpkat(*)
      DOUBLE PRECISION vargiv, epsin, epsout, varopt
      DOUBLE PRECISION X(*), soreig(*), vpmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, piw, pdw, peig
      DOUBLE PRECISION trace, som, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c     external subroutines
      EXTERNAL schurso, greps
c     
c     intrinsic functions
      INTRINSIC ABS
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw    : pointer for workspaces SCHURSO needs  (12*n)
c     pinext = piw + 12*n
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( 12*n )
c
c     pointers for double precision work space  :dwork
c     ------------------------------------------------
      peig = 1
c     peig   : pointer for vector(n), so (n) more
      pdw = peig + n
c     pdw    : pointer for workspaces SCHURSO needs  ( n*(2*n+26) )
c
c     Total size of dwork array = ( n*(2*n + 27) )
c
c-----------------------------------------------------------------------
c
c     test if vargiv in [0, 1]
      IF ( (vargiv .LT. 0.0) .OR. (vargiv .GT. 1.0) ) THEN
         info = -2003
         RETURN
      ENDIF
c
c     Schur factorization
      CALL schurso ( n, X, iwork(piw), dwork(pdw),
     &               dwork(peig), soreig, vpmat, info )
      IF (info .LT. 0) RETURN
c
c     Trace(X) 
      trace = 0.0
      DO i = 1,n
         trace = trace + soreig(i)
      ENDDO
c
c     test if |Trace(X)| < EPS
      IF ( ABS(trace) .LT. EPS ) THEN
         info = -2004
         RETURN
      ENDIF
c
c     loop until (nfact > n)
      som    = 0.0
      nbfact = 0
      varopt = - 1.0
      DO WHILE ( (varopt.LT.vargiv).and.(nbfact.LT.n) )
         nbfact = nbfact + 1
         som    = som + soreig(nbfact)
         varopt = som / trace
      ENDDO
c
c     genral case     
      IF ( nbfact .LT. n ) THEN
         epsout = (soreig(nbfact)-soreig(nbfact+1))*(1.0-EPS) - EPS
      ELSE
         epsout = soreig(1) - soreig(nbfact)
      ENDIF
c
c     case Kato-sensitivity as input
      IF (indpca .EQ. 21) THEN
        epsout = epsin
      ENDIF      
c
c     computes Kato groups for the Kato sensitivity
      CALL greps ( n, soreig, epsout, nbkgrp, grpkat, info )
c
      RETURN
      END

