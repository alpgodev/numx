c=======================================================================
c
c     subroutine IMPLRET                                    
c
c     Implied returns
c
c     ratio Sharpe = (global risk premium) / (portfolio global risk)
c                  = global risk aversion parameter
c
c     excess implied return = (global risk aversion) * Cov * w
c
c            Cov : covariance matrix
c            w   : portfolio weights
c
c     implied return = excess implied return + risk-free rate
c
c-----------------------------------------------------------------------
      SUBROUTINE implret ( n, cov, rho, w, rfr, dwork, implied, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of assets                           integer
c            cov    : covariance matrix (n*n)                     double
c            rho    : assets perf.(s) (n)                         double
c            w      : initial portfolio (n)                       double
c            rfr    : risk-free rate                              double
c
c     WORKSPACE 
c            dwork   : n                                           double
c
c     OUTPUT 
c            implied : implied return (n)                          double
c            info    : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info, i
      DOUBLE PRECISION rfr, cov(*), rho(*), w(*), implied(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER  pdwpmv
      DOUBLE PRECISION delta
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      delta = 0.0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdwpmv = 1
c     pdwpmv  : pointer for vector, so ( n ) more
c
c     Total size of dwork array = n
c
c--------------------------------------------------------------------
c     ex-ante Sharpe ratio (portfolio) - global risk aversion parameter
      CALL EXASRA ( n, cov, rho, w, rfr, delta, info)
c
c     method I, implied return (portfolio) := delta*cov*w
      CALL PMV ( n, n, cov, w, dwork(pdwpmv) )
      CALL PVX ( n, dwork(pdwpmv), delta, implied )
c
      DO i = 1,n
        implied(i) = implied(i) + rfr
      ENDDO
c
      RETURN
      END
