c=======================================================================
c
c     subroutine EXARSK                                      
c
c     Ex-ante (robust) Skewness = m3/(volatility**3)
c
c-----------------------------------------------------------------------
      SUBROUTINE EXARSK ( n, p, w, x, cov, dwork, exaske, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of value(s)                             integer
c       p       : number of asset(s)                             integer
c       w       : weights (p)                                     double
c       x       : values (n*p)                                    double
c       cov     : covariance matrix (p*p)                         double
c
c     WORKSPACE 
c       dwork   : p*(5*p + 1)                                     double
c
c     OUTPUT 
c       exaske  : ex-ante skewness                                double
c       info    : diagnostic argument                            integer
c
c     CALL 
c       EXARVO  : ex-ante variance
c       MOMENT3 : 3th moment
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION exaske
      double precision w(*), x(*), cov(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdm3
      DOUBLE PRECISION m3, vol
c
c     external subroutines
      EXTERNAL MOMENT3, EXARVO         
c
c-----------------------------------------------------------------------
c
c     initialization
      exaske = 0.0
c
c     pointers for double work space : dwork
c     -------------------------------------- 
      pdm3 = 1 
c     pdrho : pointer for average return, so (p*(5*p + 1) more
c
c     Total dwork = p*(5*p + 1)
c
c-----------------------------------------------------------------------
c
c     test number of values >= 3
      IF (n .LT. 3) THEN
        info = -2
        RETURN
      ENDIF  
c
c     covariance matrix and average return 
      CALL MOMENT3 ( n, p, w, x, dwork(pdm3), m3, info)
      IF (info .LT. 0) RETURN
c
c     ex-ante variance
      CALL EXARVO (p, cov, w, vol, info)     
      IF (info .LT. 0) RETURN
c
c     ex-ante skewness
      exaske = m3 / (vol*vol*vol)
c 
      RETURN
      END
