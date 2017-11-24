c=======================================================================
c
c     subroutine GWNPREG                                  
c
c     NPREG workspace
c
c-----------------------------------------------------------------------
c     INPUT 
c            n      : number of points                           integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE gwnpreg ( n, siwork, sdwork )
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork =  6 * ( n + 201 )
      RETURN
      END
