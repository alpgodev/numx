c=======================================================================
c
c     BFGS getwork utility
c
c     subroutines "get work" for BFGS methods
c
c------------------------------------------------------------------------
c     GWBFGSBOX
c
c     Getwork for the BFGSBOX method
c------------------------------------------------------------------------
      SUBROUTINE gwbfgsbox ( n, siwork, sdwork )
c------------------------------------------------------------------------
c     INPUT 
c            n      : number of variables                        integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = ( 2*n + 1 )
      sdwork = ( n*( n + 9 )/2 )
      RETURN
      END
c
c------------------------------------------------------------------------
c     GWBFGSBOXS
c
c     Getwork for the BFGSBOXS method
c------------------------------------------------------------------------
      SUBROUTINE gwbfgsboxs ( n, siwork, sdwork )
c------------------------------------------------------------------------
c     INPUT
c            n      : number of variables                        integer
c
c     OUTPUT
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = ( 2*n + 1 )
      sdwork = ( n*( n + 11 )/2 )
      RETURN
      END
