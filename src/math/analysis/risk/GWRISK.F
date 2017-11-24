c=======================================================================
c
c
c     Util getwork function                                  
c
c     subroutines "getwork" for rne/risk ex-ante and ex-post functions
c
c-----------------------------------------------------------------------
c
c     Ex-ante getwork 
c     GWEXAMSH  : get the size of the workspace needed by  EXAMSH
c     GWEXAMVA  : get the size of the workspace needed by  EXAMVA
c     GWEXARIR  : get the size of the workspace needed by  EXARIR
c     GWEXARKU  : get the size of the workspace needed by  EXARKU
c     GWEXARSK  : get the size of the workspace needed by  EXARSK
c     GWEXATER  : get the size of the workspace needed by  EXATER
c     GWEXANSHORTFALL : get the size of the workspace needed by  EXANSHORTFALL
c     GWEXAOMEGA : get the size of the workspace needed by  EXAOMEGA.F
c
c     Utils getwork
c     GWFUTVAL        : get the size of the workspace needed by  FUTVAL
c     GWPASTVAL       : get the size of the workspace needed by  PASTVAL
c     
c     Ex-post getwork
c     GWEXPIRA  : get the size of the workspace needed by  EXPIRA
c     GWEXPTER  : get the size of the workspace needed by  EXPTER
c     GWEXPNPVAR: get the size of the workspace needed by  EXPNPVAR 
c     GWEXPNPCVAR: get the size of the workspace needed by  EXPNPCVAR 
c     GWEXPDOWNSIDE : get the size of the workspace needed by  EXPDOWNSIDE
c     GWEXPCORRCOEF : get the size of the workspace needed by  EXPCORRCOEF
c     GWROLLINGPIES : get the size of the workspace needed by  ROLLINGPIE
c
c     GWCOKURT      : get the size of the workspace needed by COKURT
c     GWCOSKEW      : get the size of the workspace needed by COSKEW
c
c=======================================================================
c
c     GWROLLINGPIES, get the size of workspace GWROLLINGPIES
c
c-----------------------------------------------------------------------
      SUBROUTINE GWROLLINGPIES( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = 2*n
      RETURN
      END
c
c=======================================================================
c
c     GWEXPCORRCOEF, get the size of workspace GWEXPCORRCOE
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPCORRCOEF( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : number of values                                     integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = 4*n + 6
      RETURN
      END
c
c=======================================================================
c
c     GWEXAMSH, size of workspace EXAMSH (ex-ante modified Sharpe)
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXAMSH ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n      : number of values                                integer
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 0
      sdwork = p*( 7*p + 2*n + 3 )
      RETURN
      END
c
c=======================================================================
c
c     GWEXAMVA, get the size of workspace EXAMVA (ex-ante modified VaR) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXAMVA ( n, p, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, p, siwork, sdwork
      siwork = 0
      sdwork = p*( 6*p + n + 2 )  
      RETURN
      END
c
c=======================================================================
c
c     GWEXARIR, workspace need by EXARIR (ex-ante info. ratio) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXARIR ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of assets                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n*(2*n + 3) + 2 
      RETURN
      END
c
c=======================================================================
c
c     GWEXARSK, size of the workspace needed by  EXARSK (ex-ante Skewness) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXARSK ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p      : number of assets                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = 0
      sdwork = p*(5*p + 1)
      RETURN
      END
c
c=======================================================================
c
c     GWEXATER, size of the workspace needed by  EXATER (ex-ante Tracking Error) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXATER ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p      : number of assets                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n*(2*n + 3) + 2
      RETURN
      END
c
c=======================================================================
c
c     GWEXAOMEGA, size of the workspace needed by  EXAOMEGA (ex-ante Omega ratio) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXAOMEGA ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p      : number of assets                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER p, siwork, sdwork
      siwork = ( 13*p )
      sdwork = ( p*(5*p + 29))
      RETURN
      END
c
c=======================================================================
c
c     GWEXPIRA, size of the workspace needed by  EXPIRA (ex-post Info. ratio) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPIRA ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of dates                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                      integer
c       sdwork : size of double precision workspace             integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================
c
c     GWWEXPIRA, size of the workspace needed by  WEXPIRA (ex-post Info. ratio) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWWEXPIRA ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of dates                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                      integer
c       sdwork : size of double precision workspace             integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================
c
c     GWEXPTER, size of the workspace needed by  EXPTER (ex-post Tracking Err.) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPTER ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of dates                                 integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================
c
c     GWWEXPTER, size of the workspace needed by  EXPTER (ex-post Tracking Err.) 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWWEXPTER ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of dates                                 integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, siwork, sdwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================
c
c     GWEXPNPVAR
c
c     get the size of the workspace needed by  EXPNPVAR 
c     ex-post non parametric Value-at-Risk 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPNPVAR ( m, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       m      : number of non parametric point(s)               integer
c
c     OUTPUT 
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER m, sdwork, siwork
      siwork = 0
      sdwork = 2*m
      RETURN
      END
c
c=======================================================================
c
c     GWEXPNPCVAR
c
c     get the size of the workspace needed by  EXPNPCVAR 
c     ex-post non parametric conditional Value-at-Risk 
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPNPCVAR ( m, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       m      : number of non parametric point(s)               integer
c
c     OUTPUT 
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER m, sdwork, siwork
      siwork = 0
      sdwork = 2*m
      RETURN
      END
c
c======================================================================= 
c
c     GWEXPDOWNSIDE, size of the workspace needed by  DOWNSIDE (ex-post downside)
c
c-----------------------------------------------------------------------
      SUBROUTINE GWEXPDOWNSIDE ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, sdwork, siwork
      siwork = 0
      sdwork = n
      RETURN
      END
c======================================================================= 
c
c     GWWEXPDOWNSIDE, size of the workspace needed by  WDOWNSIDE (ex-post downside)
c
c-----------------------------------------------------------------------
      SUBROUTINE GWWEXPDOWNSIDE ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, sdwork, siwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================  
c  
c     GWFUTVAL, get the size of the workspace needed by  FUTVAL
c
c-----------------------------------------------------------------------
      SUBROUTINE GWFUTVAL ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, sdwork, siwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================   
c  
c     GWPASTVAL, get the size of the workspace needed by  PASTVAL
c
c-----------------------------------------------------------------------
      SUBROUTINE GWPASTVAL ( n, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, sdwork, siwork
      siwork = 0
      sdwork = n
      RETURN
      END
c
c=======================================================================   
c  
c     GWCOKURT, get the size of the workspace needed by COKURT
c
c-----------------------------------------------------------------------
      SUBROUTINE GWCOKURT ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER p, sdwork, siwork
      siwork = 0
      sdwork = p*(p + 1)
      RETURN
      END
c
c=======================================================================   
c  
c     GWCOSKEW, get the size of the workspace needed by COSKEW
c
c-----------------------------------------------------------------------
      SUBROUTINE GWCOSKEW ( p, siwork, sdwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p      : number of asset(s)                              integer
c
c     OUTPUT 
c       siwork : size of integer workspace                       integer
c       sdwork : size of double precision workspace              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER p, sdwork, siwork
      siwork = 0
      sdwork = p*(p + 1)
      RETURN
      END
