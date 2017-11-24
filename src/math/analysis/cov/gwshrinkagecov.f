c=======================================================================
c
c     GETWORK for SHRINKAGECOV
c
c=======================================================================
c
c     subroutine GWSHRINKAGECOV
c
c     SHRINKAGECOV workspace
c
c-----------------------------------------------------------------------
c     INPUT
c       nasset : size of portfolio                          integer
c     OUTPUT
c       sdwork : size of double precision workspace         integer
c
c-----------------------------------------------------------------------
      subroutine gwshrinkagecov(nasset, siwork, sdwork)
      implicit none
      integer nasset, sdwork, siwork
      siwork = 0
      sdwork = nasset*(2*nasset+1)
      return
      end
