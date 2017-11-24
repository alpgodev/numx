c=======================================================================
c
c     Asset Allocation getwork/workspace
c
c----------------------------------------------------------------------
c
c        GWALLOCTEST    : get the size of workspace need by ALLOTEST
c        GWALLOCTESTRFR : get the size of workspace need by ALLOTESTRFR
c        GWALLOCVOL     : get the size of workspace need by ALLOCVOL
c        GWALLOCVOLRFR  : get the size of workspace need by ALLOCVOLRFR
c        GWALLOCSR      : get the size of workspace need by ALLOCSR
c        GWALLOCSRRFR   : get the size of workspace need by ALLOCSRRFR
c        GWALLOCMV      : get the size of workspace need by ALLOCMV
c        GWALLOCIT      : get the size of workspace need by ALLOCIT
c        GWALLOCITRFR   : get the size of workspace need by ALLOCITRFR
c        GWALLOCVAR     : get the size of workspace need by ALLOCVAR
c        GWALLOCMDD     : get the size of workspace need by ALLOCMDD
c        GWALLOCMVTC    : get the size of workspace need by ALLOCMVTC
c        GWALLOCMVT     : get the size of workspace need by ALLOCMVT
c        GWALLOCVART    : get the size of workspace need by ALLOCVART
c        GWALLOCNAM     : get the size of workspace need by ALLOCNAM
c        GWALLOCNAMT    : get the size of workspace need by ALLOCNAMT
c        GWALLOCTE      : get the size of workspace need by ALLOCTE
c        GWALLOCRB      : get the size of workspace need by ALLOCRB
c        GWALLOCRBRFR   : get the size of workspace need by ALLOCRBRFR
c        GWMVE          : get the size of the workspace needed by  MVE
c        GWALLOCSKEW    : get the size of workspace need by ALLOCSKEW
c
c-----------------------------------------------------------------------
c
c     GWALLOCTEST, workspace size for ALLOCTEST
c
c-----------------------------------------------------------------------
      SUBROUTINE gwalloctest ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 4*n + 2*nin + neq + 1
      sdwork = n*(4*n + 16) + 3*nin + neq
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCTESTRFR, workspace size for ALLOCTESTRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwalloctestrfr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 6*n + 2*nin + neq + 13
      sdwork = n*(9*n+2*neq+2*nin+27)+7*nin+5*neq+41
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCVOL, workspace size for ALLOCVOL
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocvol ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + 2*nin + neq + 3
      sdwork = n*(2*n + nin + neq + 9) + 4*nin + 2*neq + 4
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCVOLRFR, workspace size ALLOCVOLRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocvolrfr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + 2*nin + neq + 7
      sdwork = n*(2*n + nin + neq + 13) + 4*nin + 2*neq + 12
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCSKEWGEN
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocskewgen ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, ndual, siwork, sdwork
      ndual = 2
      siwork = 2*ndual + 6 + 3*n + neq + 2*nin             
      sdwork =(ndual+21) + 4*n*n + 17*n + n*(neq + nin) + 
     & 2*(neq + nin) + 3*nin +3  
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCSR, workspace size for ALLOCSR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocsr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 7*n + 2*nin + neq + 7
      sdwork = n*(6*n+2*neq+3*nin+34)+9*nin+4*neq+27
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCSR, workspace size for ALLOCSRRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocsrrfr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 7*n + 2*nin + neq + 7
      sdwork = n*( 6*n + 2*neq + 3*nin + 35 ) + 4*neq + 9*nin + 27 
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCMV, workspace for ALLOCMV
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocmv ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + 2*nin + neq + 6
      sdwork = 2*n*(neq+nin+n+9)+4*neq+6*nin+2*n+15
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCMVT, workspace size for ALLOCMVT
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocmvt ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT :
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 13*n + 2*nin + neq + 6
      sdwork = n*( 13*n + 2*neq + 2*nin + 36 ) + 3*neq + 5*nin + 10
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCVAR, workspace size for ALLOCVAR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocvar ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + neq + 2*nin + 3
      sdwork = n*(2*n + neq + nin + 10) + 2*neq + 4*nin + 4
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCVART, workspace size for ALLOCVART
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocvart ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 13*n + neq + 2*nin + 6
      sdwork = n*(13*n + 2*neq + 2*nin + 38) + 3*neq + 5*nin + 11
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCIT, workspace size for ALLOCIT
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocit ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + neq + 2*nin + 6
      sdwork = n*(2*n + 2*nin + 2*neq + 18) + 4*neq + 6*nin + 15
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCMDD, workspace size for ALLOCMDD
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocmdd ( ndate, nasset, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       ndate  : number of dates                                 integer
c       nasset : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ndate, nasset, neq, nin, siwork, sdwork
      siwork = ( 14*nasset + 2*ndate + neq + 2*nin + 1 ) 
      sdwork = ( nasset*(8*nasset + ndate + neq + nin + 36)  
     &         + 4*( ndate + nin + 2*neq ))
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCTE, workspace size for ALLOCTE
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocte ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 16*n + neq + 2*nin + 1 
      sdwork = n*(7*n + neq + nin + 42) + 2*neq + 3*nin + 4
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCMVTC, workspace size for ALLOCMVTC
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocmvtc ( alpha, n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       alpha  : transaction cost parameter                       double
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, p, siwork, sdwork
      DOUBLE PRECISION alpha, eps
      PARAMETER ( eps = 1.E-8 )
c
c     case alpha = 0 (without trans. costs)      
      p = n
      sdwork = ( n*( 7*n + nin + neq + 40 ) + 4*nin + 2*neq )
      IF (alpha .GT. eps) THEN
        p = n + n
        sdwork = ( p*( 8*p + nin + neq + 46 ) + 4*nin + 2*neq )
      ENDIF
      siwork = ( 13*p + neq + 2*nin + 1 ) 
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCRB, workspace size for ALLOCRB
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocrb ( n, neq, nin, nbc, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = (2*nbc + 5*n + neq + 2*nin + 12)
      sdwork = (nbc*(nbc+21)/2+n*(8*n+21+3*neq+3*nin)+5*neq+7*nin+15)
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCRB, workspace size ALLOCRBRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocrbrfr ( n, neq, nin, nbc, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = 2*nbc + 5*n + neq + 2*nin + 11
      sdwork = nbc*(nbc+21)/2 + n*(8*n+2*neq+2*nin+26) 
     &       + 5*neq + 7*nin + 30
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWRISKBUDGETIT, workspace size for RISKBUDGETIT
c
c-----------------------------------------------------------------------
      SUBROUTINE gwriskbudgetit ( n, neq, nin, nbc, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c       nbc    : number of risk budgeting constraints (>1)       integer 
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer         
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, nbc, siwork, sdwork
      siwork = 3*nbc + 17*n + neq + 2*nin + 12 
      sdwork = nbc*(nbc+25)/2+n*(13*n+53+2*neq+2*nin)+6*neq+8*nin+14
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWMVE, workspace size for MVE
c
c-----------------------------------------------------------------------
      SUBROUTINE gwmve (n, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n       : number of asset(s)                              integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, mxopt, nfct, spadim, siwork, sdwork
      mxopt = 3 + n - n
      nfct = 3
      spadim = 1
      siwork=2*mxopt+7
      sdwork=mxopt*(mxopt+23)/2+nfct*(mxopt+4+2*spadim)+1
      RETURN
      END    
c-----------------------------------------------------------------------
c
c     GWALLOCMVRFR, workspace size for ALLOCMVRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocmvrfr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of asset(s)                              integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + neq + 2*nin + 7  
      sdwork = n*(2*n + neq + nin + 11) + 2*neq + 4*nin + 12
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCITRFR, workspace size for ALLOCITRFR
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocitrfr ( n, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, neq, nin, siwork, sdwork
      siwork = 3*n + neq + 2*nin + 7
      sdwork = n*(2*neq+2*nin+2*n+23)+6*nin+4*neq+31
      RETURN
      END
c-----------------------------------------------------------------------
c
c     GWALLOCSKEW, workspace size for ALLOCSKEW
c
c-----------------------------------------------------------------------
      SUBROUTINE gwallocskew ( n, p, neq, nin, siwork, sdwork )
c-----------------------------------------------------------------------
c     INPUT 
c       n      : number of dates                                 integer
c       p      : size of assets                                  integer
c       neq    : number of initial equality constraints          integer
c       nin    : number of initial inequality constraints        integer
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER n, p, neq, nin, siwork, sdwork
      siwork = 16*p + 4*nin + 2*neq + 31 + n - n
      sdwork = (p + 1)*(p + 1)*(p + nin + neq + 6)+2*nin+neq + 2*p + 7
     &       + (neq + 2*nin + 2*p + 7)*(neq + 2*nin + 2*p + 32)/2
     &       + ((2*neq + 2*nin + 2*p + 13)*(p + 1)*(p + 2))/2
     &       + (p + 1)*(15*p + 31)
      RETURN
      END    
