c=======================================================================
c
c     Workspace for perf. functions
c
c-----------------------------------------------------------------------
      SUBROUTINE gwbayesstein ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       p      : number of assets                           integer
c
c     OUTPUT 
c       siwork : size of integer workspace                  integer
c       sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 12*p + 3  
      sdwork = p*(19*p + 63)/2 + 2
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE gwjamesstein ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       p      : number of assets                           integer
c
c     OUTPUT
c       siwork : size of integer workspace                  integer
c       sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 12*p + 3
      sdwork = p*(19*p + 63)/2 + 2
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE gwimplret ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c            n      : number of assets                           integer
c
c     OUTPUT
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n
      RETURN
      END
