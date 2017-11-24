c=======================================================================
c
c     workspace for MULTIVOL methods
c
c------------------------------------------------------------------------
c
c     GWMULTIVOL                                              
c
c     Getwork for the MULTIVOL method
c
c------------------------------------------------------------------------
      SUBROUTINE gwmultivol ( n, neq, nin, nbc, siwork, sdwork )
c------------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = 2*nbc + 5*n + neq + 2*nin + 6
      sdwork = nbc*(nbc+21)/2 + n*(5*n+neq+nin+12) + 2*neq + 4*nin
      RETURN
      END
c------------------------------------------------------------------------
c
c     GWMULTIVOLRFR, getwork for the MULTIVOLRFR method
c
c------------------------------------------------------------------------
      SUBROUTINE gwmultivolrfr ( n, neq, nin, nbc, siwork, sdwork )
c------------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = 2*nbc + 5*n + neq + 2*nin + 11
      sdwork = nbc*(nbc+21)/2 + n*(8*n+2*neq+2*nin+26) 
     &       + 5*neq + 7*nin + 30
      RETURN
      END
c------------------------------------------------------------------------
c
c     GWMULTIVOLX                                             
c
c     Getwork for the MULTIVOLX method
c
c------------------------------------------------------------------------
      SUBROUTINE gwmultivolx ( n, neq, nin, nbc, siwork, sdwork )
c------------------------------------------------------------------------
c
c     INPUT 
c       n      : portfolio size                                  integer
c       neq    : number equality constraints                     integer
c       nin    : number inequality constraints                   integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = 2*nbc + 5*n + neq + 2*nin + 6
      sdwork = nbc*(nbc+19)/2 + n*(5*n+neq+nin+12) + 2*neq + 4*nin  
      RETURN
      END
