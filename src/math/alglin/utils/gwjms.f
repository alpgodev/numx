c=======================================================================
c
c     Workspace for JMS
c
c=======================================================================
c
c     subroutine GWJMS
c
c     get the size of workspace needed by JMS
c
c=======================================================================
c     INPUT :
c            n : dimension of the matrix                         integer
c
c     OUTPUT :
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
      subroutine gwjms ( n, sdwork )
      implicit none
      integer n, sdwork
      sdwork = n*n
      return
      end


