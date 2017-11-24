c=======================================================================
c
c     subroutine KATEPS
c
c     This function computes the Kato decomposition of a symmetric matrix
c     and the maximum Kato sensibility to separate choosing the next 
c     eigenvalue
c
c-----------------------------------------------------------------------
      SUBROUTINE kateps ( n, X, nsep, epsin, indpca,
     &                    iwork, dwork,
     &                    epskat, soreig, vpmat, nbkgrp, grpkat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : matrix (n*n)                                double
c            nsep   : step to separate (=1 the next)             integer
c            epsin  : Kato-sensitivity                            double
c            indpca : indicator                                  integer
c
c     WORKSPACE 
c            iwork  : vector ( 12*n )                            integer
c            dwork  : matrix ( n*(2*n+27) )                       double
c
c     OUTPUT 
c            epskat : maximum Kato sensibility                    double
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
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, nsep, nbkgrp, indpca, info
      INTEGER grpkat(*)
      DOUBLE PRECISION epsin, epskat
      DOUBLE PRECISION soreig(*), X(*), vpmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pdw, peig, i
      DOUBLE PRECISION myzero
      PARAMETER ( myzero = 1.E-30 )
c
c     external subroutines
      EXTERNAL schurso, greps
c     
c     intrinsic functions
      INTRINSIC min
c
c-----------------------------------------------------------------------
c
c     initialisations
      info = 0 
c
c     test if nsep in [1,...,n]
      IF ( ( nsep.lt.1 ).or.( nsep.gt.n ) ) THEN
         info = -2002 
         RETURN
      ENDIF
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw    : pointer for workspaces SCHURSO needs (12*n)
c     pinext = piw + ( 12*n )
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array ( 12*n )
c
c     pointers for double precision work space  : work(nasmax,*)
c     ----------------------------------------------------------
      peig = 1
c     peig   : pointer for vector(n), so (n) more
      pdw = peig + ( n )
c     pdw    : pointer for workspaces SCHURSO needs ( n*(2*n+26) )
c
c     Total size of dwork array = ( n*(2*n + 27) )
c
c-----------------------------------------------------------------------
c
c     Schur factorization
      CALL schurso ( n, X, iwork(piw), dwork(pdw),
     &               dwork(peig), soreig, vpmat, info )
      IF (info .LT. 0) RETURN
c
c     Kato-sensitivity
      IF ( nsep .LT. n ) THEN
         epskat = (soreig(nsep)-soreig(nsep+1))*(1-myzero) - myzero
      ELSE
         epskat = soreig(1)-soreig(n)
         DO i = 1,n-1
            epskat = MIN(epskat, soreig(i)-soreig(i+1))
         ENDDO
         epskat = epskat*(1 - myzero) - myzero
      ENDIF
c
c     case Kato-sensitivity as input
      IF (indpca .EQ. 11) THEN
        epskat = epsin
      ENDIF        
c
c     Kato-groups
      CALL greps ( n, soreig, epskat, nbkgrp, grpkat, info )
c
      RETURN
      END

