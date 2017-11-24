c=======================================================================
c
c     Subroutines "getwork" for Monte-Carlo methods
c
c----------------------------------------------------------------------
c
c        GWMCBASE  : get the size of the workspace needed by  MCBASE
c        GWMCBCOR  : get the size of the workspace needed by  MCBCOR
c        GWMCBHIS  : get the size of the workspace needed by  MCBHIS
c        GWMCFACT  : get the size of the workspace needed by  MCFACT
c        GWMCBASEM : get the size of the workspace needed by  MCBASEM
c        GWMCBCORM : get the size of the workspace needed by  MCBCORM
c        GWMCBHISM : get the size of the workspace needed by  MCBHISM
c        GWMCFACTM : get the size of the workspace needed by  MCFACTM
c        GWPEGBMX  : get the size of the workspace needed by  PEGBMX
c
c=======================================================================
c
c     subroutine GWMCBASE                                   
c
c     get the size of the workspace needed by  MCBASE
c
c-----------------------------------------------------------------------
c     INPUT :
c            n      : problem size                               integer
c            nstep  : number of steps to maturity                integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBASE ( n, nstep, siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, nstep, siwork, sdwork
      siwork = 0
      sdwork = ( n*(nstep + 5) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCBCOR                                    
c
c     get the size of the workspace needed by  MCBCOR
c
c-----------------------------------------------------------------------
c     INPUT :
c            n      : problem size                               integer
c            nstep  : number of steps to maturity                integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBCOR ( n, nstep, siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, nstep, siwork, sdwork
      siwork = (12*n)
      sdwork = ( n*(6*n + nstep + 27) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCBHIS                                 
c
c     get the size of the workspace needed by  MCBHIS
c
c-----------------------------------------------------------------------
c     INPUT :
c            ndate  : number of date(s)                          integer
c            n      : problem size                               integer
c            nstep  : number of steps to maturity (in years)     integer
c            h      : horizon                                    integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBHIS ( ndate, n, nstep, h,
     &                      siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, n, nstep, h, siwork, sdwork
      siwork = (12*n)
      sdwork = ( n*( 5*n + 2*ndate + 2*nstep - 2*h + 31 ) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCFACT                                    
c
c     get the size of the workspace needed by  MCFACT
c
c-----------------------------------------------------------------------
c     INPUT :
c            ndate  : number of dates                            integer
c            nasset : number of assets                           integer
c            nfact  : number of factors                          integer
c            nstep  : number of steps to maturity (in years)     integer
c            horiz  : horizon return in number of dates          integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCFACT ( ndate, nasset, nfact, nstep, horiz,
     &                      siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, nfact, nstep, horiz, siwork, sdwork
      siwork = 13*(nfact+1)
      sdwork = (
     &       nfact*( 2*ndate +3*nasset +7*nfact +2*nstep +43 -horiz )
     &     + nasset*( 5*nasset + 2*ndate + 2*nstep + 32 - 2*horiz )
     &     + ndate + 35
     &         )
      RETURN
      END
c
c=======================================================================
c
c     GWMCBASEM                                   
c
c     get the size of the workspace needed by  MCBASEM
c
c-----------------------------------------------------------------------
c     INPUT :
c            n      : size of problem                            integer
c            nstep  : number of steps to maturity                integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBASEM ( n, nstep, siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, nstep, siwork, sdwork
      siwork = 0
      sdwork = ( n*(nstep + 5) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCBCORM                                  
c
c     get the size of the workspace needed by  MCBCORM
c
c-----------------------------------------------------------------------
c     INPUT :
c            n      : size of problem                            integer
c            nstep  : number of steps to maturity (in years)     integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBCORM ( n, nstep, siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, nstep, siwork, sdwork
      siwork = (12*n)
      sdwork = ( n*(6*n + nstep + 27) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCBHISM                                   
c
c     get the size of the workspace needed by  MCBHISM
c
c-----------------------------------------------------------------------
c     INPUT :
c            ndate  : number of date(s)                          integer
c            n      : problem size                               integer
c            nstep  : number of steps to maturity                integer
c            h      : horizon                                    integer 
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCBHISM ( ndate, n, nstep, h,
     &                       siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, n, nstep, h, siwork, sdwork
      siwork = (12*n)
      sdwork = ( n*( 5*n + 2*ndate + 2*nstep - 2*h + 31 ) )
      RETURN
      END
c
c=======================================================================
c
c     GWMCFACTM                                   
c
c     get the size of the workspace needed by  MCFACTM
c
c-----------------------------------------------------------------------
c     INPUT :
c            ndate  : number of dates                            integer
c            nasset : number of assets                           integer
c            nfact  : number of factors                          integer
c            nstep  : number of steps to maturity (in years)     integer
c            horiz  : horizon return in number of dates          integer
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE GWMCFACTM ( ndate, nasset, nfact, nstep, horiz,
     &                       siwork, sdwork )
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, nfact, nstep, horiz, siwork, sdwork
      siwork = 13*(nfact+1)
      sdwork = (
     &       nfact*( 2*ndate +3*nasset +7*nfact +2*nstep +43 -horiz )
     &     + nasset*( 5*nasset + 2*ndate + 2*nstep + 32 - 2*horiz )
     &     + ndate + 35
     &         )
      RETURN
      END
c
c=======================================================================
c
c     GWPEGBMX
c
c     get the size of the workspace needed by  PEGBMX
c
c-----------------------------------------------------------------------
      SUBROUTINE GWPEGBMX ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            p      : number of assets                           integer
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 12*p
      sdwork = p*(5*p + 31)
      RETURN
      END
