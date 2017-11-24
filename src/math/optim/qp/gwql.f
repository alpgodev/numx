c=======================================================================
c
c     QL getwork utility                               
c
c     QL workspace
c
c------------------------------------------------------------------------
      SUBROUTINE gwql ( n, nceg, ncineg, siwork, sdwork )
c------------------------------------------------------------------------
c     INPUT 
c            n      : problem size                               integer
c            nceg   : number of equality constraints             integer
c            ncineg : number of inequality constraints           integer
c
c     OUTPUT 
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, nceg, ncineg, siwork, sdwork
      siwork = n
      sdwork = (3*n*n)/2 + 10*n + 2*(nceg + ncineg + 1) 
     &         + 2*n*(nceg + ncineg)
      RETURN
      END
