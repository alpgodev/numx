c=======================================================================
c
c     subroutine  GWEIGENVD   
c
c     EIGENVD workspace sizes
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEIGENVD ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : matrix size                                     integer
c
c     OUTPUT 
c       siwork : integer workspace                               integer
c       sdwork : double precision workspace                      integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      sdwork = n*n + 10*n
      siwork = 0
      RETURN
      END
