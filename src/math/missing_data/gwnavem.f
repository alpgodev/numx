c=======================================================================
c
c     subroutine GWNAVEM                                  
c
c     NAVEM workspace
c
c-----------------------------------------------------------------------
      SUBROUTINE gwnavem ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of values                                integer
c       p      : number of assets                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 2*p
      sdwork = p*(13*p + 5*n + 10)
      RETURN
      END
