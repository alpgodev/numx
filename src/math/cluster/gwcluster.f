c=======================================================================
c
c     Workspace for Cluster methods
c
c=======================================================================
c
c     subroutine GWCLUSTER
c
c-----------------------------------------------------------------------
      SUBROUTINE GWCLUSTER ( n, d, k, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT      
c       n  : number of assets                                    integer 
c       d  : dimension                                           integer
c       k  : number of cluster                                   integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, d, k, siwork, sdwork
      siwork = d
      sdwork = d*(3*d + k + 9) + n    
      RETURN
      END
