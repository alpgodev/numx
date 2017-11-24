c=======================================================================
c
c     NLS/RNLS Utilities
c
c     subroutines getwork for NLS/RNLS methods
c
c-----------------------------------------------------------------------
c
c        GWNLS    : get the size of the workspace needed by  NLS
c        GWNLSX   : get the size of the workspace needed by  NLSX
c        GWRNLS   : get the size of the workspace needed by  RNLS
c        GWRNLSX  : get the size of the workspace needed by  RNLSX
c        GWNLSB   : get the size of the workspace needed by  NLSB
c        GWRNLSB  : get the size of the workspace needed by  RNLSB
c        GWNLSS   : get the size of the workspace needed by  NLSS
c        GWRNLSS  : get the size of the workspace needed by  RNLSS
c        GWNLSBS  : get the size of the workspace needed by  NLSBS
c        GWRNLSBS : get the size of the workspace needed by  RNLSBS
c        GWNLSSX  : get the size of the workspace needed by  NLSSX
c        GWRNLSSX : get the size of the workspace needed by  RNLSSX
c        GWNLSBX  : get the size of the workspace needed by  NLSBX
c        GWRNLSBX : get the size of the workspace needed by  RNLSBX
c
c========================================================================
c
c     subroutine GWNLS
c
c     NLS workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwnls ( mxopt, nfct, spadim, liudat, ldudat,
     &                   siwork, sdwork )
c
      implicit none
      integer mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(mxopt+17)/2 + nfct*(mxopt+3+spadim) + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWNLSX
c
c     NLSX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwnlsx ( mxopt, nfct, spadim, liudat, ldudat,
     &                    siwork, sdwork )
c
      implicit none
      integer mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(mxopt+9)/2 + nfct*(mxopt+3+spadim) + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWRNLS
c
c     RNLS workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwrnls ( mxpar, nfct, spadim, liudat, ldudat,
     &                    siwork, sdwork )
c
      implicit none
      integer mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+21)/2 + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWRNLSX
c
c     RNLSX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwrnlsx ( mxpar, nfct, spadim, liudat, ldudat,
     &                     siwork, sdwork )
c
      implicit none
      integer mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(2*mxopt+19)/2 + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWNLSB
c
c     NLSB workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : integer workspace size                     integer
c            sdwork : double precision workspace size            integer
c-----------------------------------------------------------------------
      subroutine gwnlsb ( mxopt, nfct, spadim, liudat, ldudat,
     &                    siwork, sdwork )
c
      implicit none
      integer mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(mxopt+13)/2 + nfct*(mxopt+3+spadim) + ldudat )
      return
      end
c      
c=======================================================================
c
c     subroutine GWNLSBX
c
c     NLSBX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwnlsbx ( mxopt, nfct, spadim, liudat, ldudat,
     &                     siwork, sdwork )
c
      implicit none
      integer mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*( mxopt + 11 )/2 
     &         + nfct*( mxopt + 3 + spadim ) + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWRNLSB
c
c     RNLSB workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwrnlsb ( mxpar, nfct, spadim, liudat, ldudat,
     &                     siwork, sdwork )
c
      implicit none
      integer mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+19)/2 + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWNLSS
c
c     NLSS workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
c
      subroutine gwnlss ( mxopt, nfct, spadim, liudat, ldudat,
     &                    siwork, sdwork )
c
      implicit none
      integer mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(3*mxopt+21)/2 + nfct*(mxopt+3+spadim) + ldudat )
      return
      end
c
c=======================================================================
c
c     subroutine GWRNLSS
c
c     RNLSS workspace
c
c-----------------------------------------------------------------------
c     INPUT 
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE gwrnlss ( mxpar, nfct, spadim, liudat, ldudat,
     &                     siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+19)/2 + (mxpar*(mxpar+2))
     &           + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWNLSBS
c
c     NLSBS workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE gwnlsbs( mxopt, nfct, spadim, liudat, ldudat,
     &                    siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(3*mxopt+17)/2 + nfct*(mxopt+3+spadim) + ldudat )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWRNLSBS
c
c     RNLSBS workspace
c
c-----------------------------------------------------------------------
c     INPUT 
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      SUBROUTINE gwrnlsbs ( mxpar, nfct, spadim, liudat, ldudat,
     &                      siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+19)/2 + (mxpar*(mxpar+2))
     &           + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWNLSSX
c
c     NLSSX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : integer workspace size                     integer
c            sdwork : double precision workspace size            integer
c-----------------------------------------------------------------------
      SUBROUTINE gwnlssx ( mxopt, nfct, spadim, liudat, ldudat,
     &                     siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxopt, nfct, spadim, liudat, ldudat, siwork, sdwork
      siwork = ( 2*mxopt + 6 + liudat )
      sdwork = ( mxopt*(3*mxopt+13)/2 + nfct*(mxopt+3+spadim) + ldudat )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWRNLSSX
c
c     RNLSSX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : integer workspace size                     integer
c            sdwork : double precision workspace size            integer
c-----------------------------------------------------------------------
      SUBROUTINE gwrnlssx ( mxpar, nfct, spadim, liudat, ldudat,
     &                      siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+9)/2 + (mxpar*(mxpar+2))
     &           + nfct*(7*nfct+mxpar+spadim+33)
     &           + 2 + ldudat )
      RETURN
      END
c
c=======================================================================
c
c     subroutine GWRNLSBX
c
c     RNLSBX workspace
c
c-----------------------------------------------------------------------
c     INPUT :
c            mxpar  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            spadim : space dimension                            integer
c            liudat : size of integer user data for SIMEXT       integer
c            ldudat : size of double prec. user data for SIMEXT  integer
c
c     OUTPUT :
c            siwork : integer workspace size                     integer
c            sdwork : double precision workspace size            integer
c-----------------------------------------------------------------------
c
      SUBROUTINE gwrnlsbx ( mxpar, nfct, spadim, liudat, ldudat,
     &                      siwork, sdwork )
c
      IMPLICIT NONE
      INTEGER mxpar, nfct, spadim, liudat, ldudat, siwork, sdwork, mxopt
      mxopt = mxpar + 1
      siwork = ( 2*mxopt + 6 + 12*nfct + liudat )
      sdwork = ( mxopt*(mxopt+13)/2 
     &           + nfct*(7*nfct + mxpar + spadim + 33)
     &           + 2 + ldudat )
      RETURN
      END
