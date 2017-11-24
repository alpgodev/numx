c=======================================================================
c
c     QP getwork utility
c
c     QP subroutines
c
c------------------------------------------------------------------------
      SUBROUTINE GWQP ( n, nceg, ncineg, siwork, sdwork )
c------------------------------------------------------------------------
c     INPUT 
c            n      : problem size                               integer
c            nceg   : number of equality constraints             integer
c            ncineg : number of inequality constraints           integer
c
c     OUTPUT 
c            siwork : integer workspace size                     integer
c            sdwork : double precision workspace size            integer
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, nceg, ncineg, siwork, sdwork
      siwork = ( 3*n + 2*ncineg + nceg + 1 )
      sdwork = ( n*n + 6*n + 2*ncineg )
      RETURN
      END
