c=======================================================================
c
c     subroutine CHTMCM
c
c     Computes the mean of each column of a vectorized matrix
c     with a poor historic (some unknown - "unborn" in fact - values 
c     at the beginning of historic only)
c
c-----------------------------------------------------------------------
      SUBROUTINE chtmcm ( nlmat, ncmat, mat, hole, ages, vect, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows (= number of dates) 
c                     of the matrix mat                          integer
c            ncmat  : number of columns of the matrix mat 
c                     (= number of assets)                       integer
c            mat    : matrix as vector (nlmat*ncmat)              double
c            hole   : stands for unknown value                    double
c            ages   : vector (ncmat) of "ages of the assets
c                     (=nb of known values in the historic)      integer
c
c     OUTPUT 
c            vect   : mean of each column, vector(ncmat)          double
c            info   : diagnostic argument                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER nlmat, ncmat, info, ages(*)
      DOUBLE PRECISION hole, mat(nlmat,*), vect(*)
c
c     local variables
      INTEGER j, k, ind1st
      DOUBLE PRECISION sum, eps
      PARAMETER( eps = 1.E-5 )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      DO j = 1,ncmat
         IF (ages(j) .LT. 1) THEN
            info = -2
            RETURN
         ENDIF
         ind1st = nlmat-ages(j) + 1
c         
c        sum
         sum = 0.
         DO k = ind1st,nlmat
            IF (      (mat(k,j) .GE. hole - eps) 
     &          .and. (mat(k,j) .LE. hole + eps)) THEN
c
c              hole value
               info = -4005
               RETURN
            ENDIF
            sum = sum + mat(k,j)
         ENDDO
         vect(j) = sum / (ages(j))
      ENDDO
c
      RETURN
      END
