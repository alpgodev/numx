c=======================================================================
c
c     subroutine PCA                                         
c
c     PCA (Principal Component Analysis)
c
c-----------------------------------------------------------------------
      SUBROUTINE pca ( ndate, nasset, dasset, indpca,
     &                 nfask, varmin, epsin,
     &                 iwork, dwork, nfatra, varian, covcor, 
     &                 specor, eigvec, epskat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            ndate  : full estimation period                     integer
c            nasset : size of problem                            integer
c            dasset : returns  matrix(ndate*nasset)               double
c            indpca : indicator                                  integer
c                = 1   factor tracking with nfask-factor(s) 
c                  11  (epskat constant input)
c                  12  (epskat constant computed)
c                = 2 factor tracking for a level of explicative variance
c                  11  (epskat constant input)
c                  12  (epskat constant computed)
c            nfask  : asked number of tracked factors            integer
c            varmin : minimun variance to reach                   double
c            epsin  : Kato-sensitivity in input                   double
c
c     WORKSPACE 
c            iwork  : vector( 13*nasset )                        integer
c            dwork  : vector( nasset *(ndate+3*nasset+27) )     double
c
c     OUTPUT 
c            nfatra : number of tracked factors                  integer
c            epskat : computed Kato sensibility                   double
c            varian : variance                                    double
c            covcor : corrected covariance matrix(nasset*nasset)  double
c            specor : corrected spectrum vector(nasset)           double
c            eigvec : eigenvectors matrix(nasset*nasset)          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c        COVMVM    : covariance matrix (equally weighted) and mean vector
c        KATVAR  : Kato decomposition of a symmetric matrix
c                  the maximum Kato sensibility and the number factors
c                  for a given variance
c        KATEPS  : computes the Kato decomposition of a symmetric matrix
c                  and the maximum Kato sensibility
c                  to separate choosing eigenvalue of the next
c        OMCDMCT : computing M*D*M'
c                  (M square matrix n*n, D vector n of diagonal matrix,
c                   gives square matrix n*n )
c        ZMCMC   : symmetrises a vectorized full square matrix
c                  and gives a vectorized full square matrix A=(A+A')/2
c        VARFAC  : computes the variance of n first factors
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER ndate, nasset, indpca, nfask, nfatra, ngrkat, info
      DOUBLE PRECISION varmin, epsin, epskat, varian
      DOUBLE PRECISION dasset(ndate,*), covcor(*), specor(*), eigvec(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
      INTEGER iwork(*)
c
c     local parameters
      INTEGER i, ig, iev, nbev, pigk, piw, pdw, pvecw, pmatw
      DOUBLE PRECISION mean, sum
c
c-----------------------------------------------------------------------
c
c     test indpca
      IF ((indpca .NE. 11).AND.(indpca .NE. 12)
     &    .AND.(indpca .NE. 21).AND.(indpca .NE. 22)) THEN
         info = -4001
         RETURN
      ENDIF
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pigk = 1
c     pigk  : pointer for Kato groups,
      piw = pigk + nasset
c     piw   : pointer for workspaces,
c             KATVAR needs (12*nasset)
c             KATEPS needs (12*nasset)
c             so  (12*nasset) more
c
c     Total size of iwork array = (13*nasset)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pvecw = 1
c     pvecw : pointer for work vector, so (nasset) more
      pmatw = pvecw + (nasset)
c     pdata : pointer for work matrix, so (nasset*nasset) more
      pdw   = pmatw + (nasset*nasset)
c     pdw   : pointer for workspaces,
c             COVMVM   needs ( ndate*nasset )
c             KATVAR needs ( nasset*(2*nasset+27) )
c             OMCDMCT needs ( nasset*nasset )
c             so  ( nasset*(ndate+2*nasset+27) ) more
c
c     Total size of dwork array = (nasset) + (nasset*nasset)
c                               + ( nasset*(ndate+2*nasset+27)
c              ( nasset *( ndate + 3*nasset + 28 ) )
c
c-----------------------------------------------------------------------
c
c     initial covariance matrix (empiric equally weighted)
      CALL covmvm ( ndate, nasset, dasset, dwork(pdw),
     &            dwork(pvecw), covcor, info)
      IF (info .LT. 0) RETURN
c
c     computing Kato sensibily and number of tracked factors
c     to approach as best the asked variance
      IF ((indpca .EQ. 21) .AND. (indpca .EQ. 22)) THEN
         CALL katvar ( nasset, covcor, varmin, epsin, indpca,
     &                 iwork(piw), dwork(pdw), varian, nfatra,
     &                 epskat, specor, eigvec,
     &                 ngrkat, iwork(pigk), info )
         IF (info .LT. 0) RETURN
      ELSE
         nfatra = nfask
         CALL kateps ( nasset, covcor, nfatra, epsin, indpca,
     &                 iwork(piw), dwork(pdw), epskat,
     &                 specor, eigvec,
     &                 ngrkat, iwork(pigk), info )
         IF (info .LT. 0) RETURN
      ENDIF
c
c     averages eigenvalues in Kato groups
      iev = 0
      DO ig = 1,ngrkat
         nbev = iwork(pigk+ig-1)
         sum  = 0.
c
         DO i = 1,nbev
            sum = sum + specor(iev+i)
         ENDDO
c
         mean = sum / nbev
c
         DO i = 1,nbev
            specor(iev+i) = mean
         ENDDO
c
         iev = iev + nbev
      ENDDO
c
c     covariance matrix = eigvec*specor*eigvec'
      CALL OMCDMCT ( nasset, eigvec, specor, dwork(pdw),
     &               dwork(pmatw) )
c
c     construct a symmetry matrix (for computer precision)
      CALL ZMCMC ( nasset, dwork(pmatw), covcor )
c
c     variance of the factors (alglin/utils/utstat.f)
      CALL VARFAC ( nasset, specor, nfatra, varian, info )
c
      RETURN
      END

