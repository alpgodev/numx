c=======================================================================
c
c     GETWORK for cmc functions
c
c-----------------------------------------------------------------------
c
c        List of getwork (workspaces) functions:
c
c        GWCALEPSMEDIAN   : get the size of the workspaces needed by CALEPSMEDIAN
c        GWSCHUR   : get the size of the workspaces needed by SCHUR
c        GWSCHURPCA: get the size of the workspaces needed by SCHURPCA
c        GWSCHURS  : get the size of the workspaces needed by SCHURS
c        GWSCHURO  : get the size of the workspaces needed by SCHURO
c        GWSCHURSO : get the size of the workspaces needed by SCHURSO
c        GWKAT     : get the size of the workspaces needed by KAT
c        GWKATP    : get the size of the workspaces needed by KATP
c        GWKATAVE  : get the size of the workspaces needed by KATAVE
c        GWKATAVEP : get the size of the workspaces needed by KATAVEP
c        GWKATEPS  : get the size of the workspaces needed by KATEPS
c        GWKATVAR  : get the size of the workspaces needed by KATVAR
c        PROJECT   : get the size of the workspaces needed by PROJECT
c        PROJECTS  : get the size of the workspaces needed by PROJECTS
c        PROJECTA  : get the size of the workspaces needed by PROJECTA
c        GWCHOL    : get the size of the workspaces needed by CHOL
c        GWRCHO    : get the size of the workspaces needed by RCHO
c        GWEPSINIT : get the size of the workspaces needed by EPSINIT
c        GWKATOBLOCK  :
c        GWKATOBLOCKX :
c        GWCALEPSVAR  :
c
c=======================================================================
c
c     subroutine GWCALEPSMEDIAN
c
c     get the size of the workspaces needed by CORM
c
c-----------------------------------------------------------------------
      SUBROUTINE GWCALEPSMEDIAN ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n  : estimation period                                   integer
c
c     OUTPUT 
c       siwork : integer workspace                               integer
c       sdwork : double precision workspace                      integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, sdwork, siwork
      sdwork = 2*n*n + 9*n
      siwork = 0
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSCHUR
c
c     get the size of the workspaces needed by SCHUR
c
c-----------------------------------------------------------------------
      SUBROUTINE GWSCHUR ( n, siwork, sdwork )
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
      siwork = n
      sdwork = n*(n + 6)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSCHURPCA
c
c     get the size of the workspaces needed by SCHURPCA
c
c-----------------------------------------------------------------------
      SUBROUTINE GWSCHURPCA ( n, siwork, sdwork )
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
      siwork = n
      sdwork = 2*n*(n + 7)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSCHURS
c
c     get the size of the workspaces needed by SCHURS
c
c-----------------------------------------------------------------------
      SUBROUTINE GWSCHURS ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSCHURO
c
c     get the size of the workspaces needed by SCHURO
c
c-----------------------------------------------------------------------
      SUBROUTINE GWSCHURO ( n, siwork, sdwork )
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
      siwork = n
      sdwork = n*(2*n + 6)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSCHURSO
c
c     get the size of the workspaces needed by SCHURSO
c
c-----------------------------------------------------------------------
      SUBROUTINE GWSCHURSO ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of matrix                                  integer
c
c     OUTPUT 
c       siwork : integer workspace                               integer
c       sdwork : double precision workspace                      integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 12*n
      sdwork = n*(2*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKAT
c
c     get the size of the workspaces needed by KAT
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKAT ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(2*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATP
c
c     get the size of the workspaces needed by KATP
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKATP ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of matrix                                  integer
c
c     OUTPUT 
c       siwork : integer workspace                               integer
c       sdwork : double precision workspace                      integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 12*n
      sdwork = n*(2*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATAVE
c
c     get the size of the workspaces needed by KATAVE
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKATAVE ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(4*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATAVEP
c
c     get the size of the workspaces needed by KATAVEP
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKATAVEP ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(4*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATEPS
c
c     get the size of the workspaces needed by KATEPS
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKATEPS ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(2*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATVAR
c
c     get the size of the workspaces needed by KATVAR
c
c-----------------------------------------------------------------------
      SUBROUTINE GWKATVAR ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(2*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWPROJECT
c
c     get the size of the workspaces needed by PROJECT
c
c-----------------------------------------------------------------------
      SUBROUTINE GWPROJECT ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(3*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWPROJECTS
c
c     get the size of the workspaces needed by PROJECTS
c
c-----------------------------------------------------------------------
      SUBROUTINE GWPROJECTS ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(3*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWPROJECTA
c
c     get the size of the workspaces needed by PROJECTA
c
c-----------------------------------------------------------------------
      SUBROUTINE GWPROJECTA ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(3*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCHOL
c
c     get the size of the workspaces needed by CHOL
c
c-----------------------------------------------------------------------
      SUBROUTINE GWCHOL ( n, siwork, sdwork )
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
      siwork = 0
      sdwork = n*n
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWRCHO
c
c     get the size of the workspaces needed by RCHO
c
c-----------------------------------------------------------------------
      SUBROUTINE GWRCHO ( n, siwork, sdwork )
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
      siwork = 12*n
      sdwork = n*(4*n + 27)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWEPSINIT
c
c     get the size of the workspaces needed by EPSINIT
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEPSINIT ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of eigenvalues                           integer
c
c     OUTPUT 
c       siwork : integer workspace                               integer
c       sdwork : double precision workspace                      integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 12*n
      sdwork = n*(2*n + 26)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATOBLOCK
c
c-----------------------------------------------------------------------
      SUBROUTINE gwkatoblock ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : matrix size                                     integer
c
c     OUTPUT
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 14*n + 3
      sdwork = 5*(n*(n + 1)/2) + n*(8*n + 30)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWKATOBLOCKX
c
c-----------------------------------------------------------------------
      SUBROUTINE gwkatoblockx ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : matrix size                                     integer
c
c     OUTPUT
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 14*n + 3
      sdwork = 5*(n*(n + 1)/2) + n*(8*n + 30)
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWCALEPSVAR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcalepsvar ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : number of eigenvalues                           integer
c
c     OUTPUT
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = n
      sdwork = n*(2*n + 17)
      RETURN
      END
