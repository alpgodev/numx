c=======================================================================
c
c     subroutine KATAVE
c
c     This function calibrates a symmetric matrix by average Kato
c
c-----------------------------------------------------------------------
      SUBROUTINE katave ( n, X, epskat, iwork, dwork,
     &                    Y, spectr, eigsor, eigave, 
     &                    nbkgrp, grpkat, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : matrix size                                integer
c            X      : symmetric matrix (n*n)                      double
c            epskat : Kato sensibility parameter                  double
c
c     WORKSPACE 
c            iwork  : vector ( 12*n )                            integer
c            dwork  : vector ( n*(4*n + 26) )                     double
c
c     OUTPUT 
c            Y      : calibrate matrix (n*n)                      double
c            spectr : initial spectrum (n)                        double
c            eigsor : sorted eigenvalues (n)                      double
c            eigave : average eigenvalues (n)                     double
c            nbkgrp : number of Kato groups                      integer
c            grpkat : Kato groups (n) output size (nbkgrp)       integer
c            info   : diagnostic argument                        integer
c
c     CALL   
c            KAT     : Kato decomposition of a sym. full vectorized matrix
c            OMCDMCT : computing M*D*M'
c                     (M square matrix n*n, D vector n of diagonal matrix,
c                      gives square matrix n*n )
c            ZMCMC   : symmetrises a vectorized full square matrix
c                      and gives a vectorized full square matrix A=(A+A')/2
c
c     METHOD 
c            for Kato groups vector    
c            for an example of eigenvalues vector ( size = 7 ) :   
c            eigval = [ 5.2, 5.1, 3.5, 3.4, 3.3, 1.7 1.3 ]    
c            and eps = 0.15, the routine computes    
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
      DOUBLE PRECISION X(n,*), Y(n,*), spectr(*), eigsor(*), eigave(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER iev, nbev, ig, i, piw, pmat, pvs, pdw
      DOUBLE PRECISION som, mean
c
c     external subroutines
      EXTERNAL kat, OMCDMCT, ZMCMC
c
c-----------------------------------------------------------------------
c
c     initialization 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for KAT who needs  ( 12*n )
c     pinext = piw + 12*n
c     pinext : pointer for the next iwork array
c
c     Total size of iwork array = ( 12*n )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pmat = 1
c     pmat  : pointer for a work matrix , so (n*n) more
      pvs = pmat + ( n * n )
c     pvs   : pointer eigenvectors matrix, so (n*n) more
      pdw = pvs + ( n * n )
c     pdw   : pointer for workspaces,
c             KAT needs ( n*(2*n + 26) )
c             OMCDMCT needs (n*n)
c             so ncov*(2*n+26) more
c     pdnext = pdw + n*(2*n+26)
c     pdnext: pointer of the next dwork array
c
c     Total size of dwork array = ( n*(4*n+26) )
c
c-----------------------------------------------------------------------
c
c     Kato decomposition
      CALL kat ( n, X, epskat, iwork(piw), dwork(pdw),
     &           spectr, eigsor, dwork(pvs), nbkgrp, grpkat, info )
      IF (info .LT. 0) RETURN
c
c     averages eigenvalues in Kato groups
      iev = 0
      DO ig = 1,nbkgrp
         nbev = grpkat(ig)
         som  = 0.0
         DO i = 1,nbev
            som = som + eigsor(iev+i)
         ENDDO
         mean = som / nbev
         DO i = 1,nbev
            eigave(iev+i) = mean
         ENDDO
         iev = iev + nbev
      ENDDO
c
c     reconstruction of the matrix = Z*T*Z'
      CALL OMCDMCT ( n, dwork(pvs), eigave, dwork(pdw), dwork(pmat) )
c
c     assures the symetry of the matrix (computer precision)
      CALL ZMCMC ( n, dwork(pmat), Y )
c
      RETURN
      END

