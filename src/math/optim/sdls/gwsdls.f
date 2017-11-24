c=======================================================================
c
c     SDLS getwork utilities                               
c
c-----------------------------------------------------------------------
c
c        GWSDLS      : SDLS workspace
c        GWSDLSCOR   : SDLSCOR workspace
c        GWSDLSCE    : SDLSCE workspace
c        GWSDLSG     : SDLSG workspace
c        GWSDLSGEN   : SDLSG workspace
c        GWSDLSCI    : SDLSCI workspace
c        GWSDLSTRACE : SDLSTRACE workspace
c
c=======================================================================
c
c     subroutine GWSDLS
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdls ( n, mct, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : dimension of the matrix                    integer
c            mct    : number of constraints                      integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, mct, siwork, sdwork
      siwork = ( 2*mct + 12*n + 3 )
      sdwork = ( mct*(mct+25)/2 + (2*mct+5)*(n*(n + 1)/2)
     &                          + (n*(4*n + 27)) )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSCOR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlscor ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : matrix size                                integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = ( 14*n + 3 )
      sdwork = ( n*(7*n + 40) )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSCE
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlsce ( n, imatce, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : dimension of the matrix                    integer
c            imatce : symmetric matrix (n*n)                     integer
c                     equal constraints indicators
c                     only the lower part is used
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
c     i/o arguments
      INTEGER n, imatce(n,*), siwork, sdwork
c     local parameters
      INTEGER mcte, i, j
c     computing the number of equal constraints
      mcte = 0
      DO i = 1,n
         DO j = 1,i
            IF ( imatce(i,j) .NE. 0 ) mcte = mcte + 1
         ENDDO
      ENDDO
      siwork = ( 2*mcte + 12*n + 3 )
      sdwork = (   mcte*(n*n)
     &           + mcte*(mcte + 27)/2
     &           + (2*mcte + 5)*(n*(n + 1)/2)
     &           + (n*(4*n + 27)) )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSG
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlsg ( n, mcte, mcti, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : matrix size                                integer
c            mcte   : number of equal constraints                integer
c            mcti   : number of inequal constraints              integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, mcte, mcti, siwork, sdwork
      siwork = ( 2*mcte + 4*mcti + 12*n + 5 )
      sdwork = (   (mcte + 2*mcti)*(mcte + 2*mcti + 25)/2
     &           + (2*mcte + 2*mcti + 5)*(n*(n + 1)/2)
     &           + (n*(4*n + 27))  )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSGEN
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlsgen ( n, mcte, mcti, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : matrix size                                integer
c            mcte   : number of equal constraints                integer
c            mcti   : number of inequal constraints              integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, mcte, mcti, siwork, sdwork
      CALL gwsdlsg ( n, mcte, mcti, siwork, sdwork )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSCI
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlsci ( n, imatc, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : matrix size                                integer
c            imatc  : symmetric matrix (n*n)                     integer
c                     constraints indicators (lower part is used)
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
c     i/o arguments
      INTEGER n, siwork, sdwork, imatc(n,*)
c     local variables
      INTEGER mcte, mcti, i, j
c     number of equal and inequal constraints
      mcte = 0
      mcti = 0
      DO i = 1,n
         DO j = 1,i
            IF ( imatc(i,j) .EQ. 1 ) mcte = mcte + 1
            IF ( imatc(i,j) .EQ. 2 ) mcte = mcte + 1
            IF ( imatc(i,j) .EQ. 3 ) mcti = mcti + 1
            IF ( imatc(i,j) .EQ. 4 ) mcti = mcti + 1
         ENDDO
      ENDDO
      siwork = ( 2*mcte + 4*mcti + 12*n + 5 )
      sdwork = (  (mcti+mcte)*(n*n)
     &          + (mcti) + (mcti) + mcti*(n*n)
     &          + (mcte+2*mcti)*((mcte+2*mcti+27)/2)
     &          + (2*mcte+2*mcti+5)*(n*(n+1)/2)
     &          + (n*(4*n + 27))  )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWSDLSTRACE
c
c-----------------------------------------------------------------------
      SUBROUTINE gwsdlstrace ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            n      : dimension of the matrix                    integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 12*n + 5
      sdwork = n*(7*n + 27) + 6 
      RETURN
      END
