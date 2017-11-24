c=======================================================================
c
c     subroutine KAT
c
c     Kato decomposition of a symmetric matrix
c
c----------------------------------------------------------------------
      SUBROUTINE kat ( n, X, epskat, iwork, dwork,
     &                 eigval, soreig, vpmat, nbkgrp, grpkat, info )
c----------------------------------------------------------------------
c
c     INPUT 
c            n      : size of matrix X                           integer
c            X      : matrix to decompose (n*n)                   double
c            epskat : Kato sensibility parameter                  double
c
c     WORKSPACE 
c            iwork  : vector ( 12*n )                            integer
c            dwork  : matrix ( n*(2*n+26) )                       double
c
c     OUTPUT 
c            eigval : eigenvalues (n)                             double
c            soreig : sorted eigenvalues (n)                      double
c            vpmat  : ordered eigenvectors (n*n)                  double
c            nbkgrp : number of Kato groups                      integer
c            grpkat : Kato groups (n) output size (nbkgrp)       integer
c            info   : diagnostic argument                        integer
c
c     CALL   
c            SCHURSO: SCHUR decomposition of a symmetric matrix
c                     sorts eigenvalues and eigenvectors (decrease order)
c             GREPS : computes groups for a sensibility parameter
c
c     METHOD 
c            for Kato groups vector    
c            for an example of eigenvalues vector ( size = 7 ) :   
c            eigval = [ 5.2, 5.1, 3.5, 3.4, 3.3, 1.7 1.3 ]    
c            and epskat = 0.15, the routine computes    
c            nbkgrp = 4    
c            grpvec = [ 2, 3, 1, 1 ]    
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, nbkgrp, info
      INTEGER grpkat(*)
      DOUBLE PRECISION epskat
      DOUBLE PRECISION eigval(*), soreig(*), x(*), vpmat(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     external subroutines
      EXTERNAL schurso, greps
c
c     SCHUR decomposition of a symmetric matrix
c     sort eigenvalues and eigenvectors (decrease order)
      CALL schurso ( n, X, iwork, dwork, eigval, soreig, vpmat, info )
      IF (info .lt. 0) RETURN
c
c     Kato epskat-groups construction
      CALL greps ( n, soreig, epskat, nbkgrp, grpkat, info )
c
      RETURN
      END

