c=======================================================================
c
c     Workspace for APT methods
c
c     GWCALAPT : CALAPT workspace
c     GWAPTVOL : APTVOL workspace
c     GWAPTCST : APTCST workspace
c
c=======================================================================
c
c     GWCALAPT, get the size of the workspace needed by  CALAPT
c
c-----------------------------------------------------------------------
      SUBROUTINE gwcalapt ( ndate, nasset, nfact, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : number of dates                            integer
c            nasset : size of assets                             integer
c            nfact  : number of factors                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, nfact, siwork, sdwork
      siwork = nfact+1
      sdwork = (nfact + 1)*( ndate + 2*nasset + nfact + 2 )
      RETURN
      END
c=======================================================================
c
c     GWAPTVOL, get the size of the workspace needed by  APTVOL
c
c-----------------------------------------------------------------------
      SUBROUTINE gwaptvol ( ndate, nasset, nfact, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : number of dates                            integer
c            nasset : size of assets                             integer
c            nfact  : number of factors                          integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, nfact, siwork, sdwork
      siwork = nfact+1
      sdwork = (nfact+1)*( ndate + 2*nasset + nfact + 2 ) + ndate*nasset
      RETURN
      END
c=======================================================================
c
c     GWAPTCST, get the size of the workspace needed by  APTCST
c
c-----------------------------------------------------------------------
      SUBROUTINE gwaptcst ( ndate, nfact, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c            ndate  : number of dates                            integer
c            nfact  : number of factors                          integer
c            neq    : number of equality constraints             integer
c            nin    : number of inequality constraints           integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nfact, neq, nin, siwork, sdwork, p
      p = nfact + 1
      siwork = 3*p + 2*nin + neq + 1
      sdwork = p*(2*p + ndate + 9) + 3*nin + neq
      RETURN
      END
