c=======================================================================
c
c     Workspace for cov
c
c-----------------------------------------------------------------------
c
c        List of getwork (workspaces) functions:
c
c        GWCORM          : CORM workspace
c        GWCOVAXP        : COVEXP workspace
c        GWCORVOLM       : COVM workspace
c        GWCOVM          : COVM workspace
c        GWCOVL          : COVL workspace
c        GWCOVLACK       : COVLACK workspace
c        GWCORLACK       : CORLACK workspace
c        GWCOVFILTERING  : COVFILTERING workspace
c        GWCOVFILTERING3 : COVFILTERING3 workspace
c
c=======================================================================
c
c     subroutine GWCORM
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcorm ( ndate, nasset, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : estimation period                          integer
c            nasset : size of portfolio                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, siwork, sdwork
      siwork = 0
      sdwork = ( nasset*(ndate+nasset+2) )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCORVOLM
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcorvolm GWCORVOLM ( ndate, nasset, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : estimation period                          integer
c            nasset : size of portfolio                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, siwork, sdwork
      siwork = 0
      sdwork = ( nasset*(ndate+nasset) )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCOVM
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovm ( ndate, nasset, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : estimation period                          integer
c            nasset : size of portfolio                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, siwork, sdwork
      siwork = 0
      sdwork = ( ndate*nasset + nasset)
      RETURN
      END
c      
c=======================================================================
c
c     subroutine GWCOVEXP
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovexp ( ndate, nasset, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : estimation period                          integer
c            nasset : size of portfolio                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, siwork, sdwork
      siwork = 0
      sdwork = ( (ndate + 1)*nasset )
      RETURN
      END   
c
c=======================================================================
c
c     subroutine GWCOVL
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovl ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : number of values                           integer
c            p      : number of assets                           integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 0
      sdwork = ( p*(1 + n + p) + 2*n )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCOVLACK
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovlack ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of date(s)                               integer
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 0
      sdwork = ( n*(p + 4) + 6 )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCORLACK
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcorlack ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of date(s)                               integer
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 0
      sdwork = p*(n + p + 1) + 4*n + 6
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCOVFILTERING
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovfiltering ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       p      : matrix size                                     integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 12*p + 3
      sdwork = (p*(15*p + 61))/2 + 1
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCOVFILTERING3
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcovfiltering3 ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       p      : matrix size                                     integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 12*p + 3
      sdwork = (p*(15*p + 61))/2 + 1
      RETURN
      END
