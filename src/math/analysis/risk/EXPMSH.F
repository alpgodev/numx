c=======================================================================
c
c     subroutine EXPMSH                                      
c
c     Ex-post Cornish-Fisher (or modified) Sharpe Ratio
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPMSH( n, x, rfr, beta, expms, info)
c-----------------------------------------------------------------------
c
c     INPUT :
c         n      : number of return(s)                           integer
c         x      : return(s) (n)                                  double
c         beta   : quantile                                       double
c         rfr    : risk-free rate                                 double
c       
c     OUTPUT :
c         expms  : ex-post modified Sharpe ratio                  double
c         info   : error code
c
c     CALL :  
c        EXPRET  : ex-post mean return
c        EXPMVA  : ex-post modified Value-at-Risk
c        MODSHA  : modified Sharpe ratio
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o   
      INTEGER n, info   
      double precision rfr, beta, expms      
      DOUBLE PRECISION x(*)
c
c     local variables
      DOUBLE PRECISION expmv, mexpr      
c
c     external subroutines
      EXTERNAL EXPRET, EXPMVA, MODSHA
c
c----------------------------------------------------------------------
c
c     initialization
      expms = 0.0
      info  = 0
c
c     ex-post mean return 
      CALL EXPRET ( n, x, mexpr, info )
      IF (info .lt. 0) RETURN
c
c     ex-post modified Value-at-Risk
      CALL EXPMVA( n, x, beta, expmv, info )
      IF (info .lt. 0) RETURN
c
c     ex-post modified Sharpe ratio
      CALL MODSHA(mexpr, rfr, expmv, expms, info)
      IF (info .lt. 0) RETURN
c      
      RETURN
      END
     
