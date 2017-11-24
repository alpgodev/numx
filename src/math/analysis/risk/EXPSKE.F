c=======================================================================
c
c     subroutine EXPSKE                                      
c
c     Ex-post Skewness
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPSKE ( n, x, expsk, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c          n       : size of x (n > 1)                           integer
c          x       : returns (n)                                  double
c
c     OUTPUT 
c          expsk   : ex-post Skewness                             double
c          info    : diagnostic argument                         integer
c
c     CALL   
c          EXPRET  : ex-post return
c          EXPVAR  : ex-post robust variance
c          SKEWN   : ex-post skewness
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION expsk
      DOUBLE PRECISION x(*)
c
c     local variables      
      DOUBLE PRECISION expr, var
c
c     external subroutines
      EXTERNAL EXPRET, EXPVAR, SKEWN      
c
c-----------------------------------------------------------------------
c
c     initialization
      expsk = 0.0
      info  = 0
c
c     ex-post return
      CALL EXPRET ( n, x, expr, info )
      IF (info .lt. 0) RETURN
c
c     ex-post variance
      CALL EXPVAR ( n, x, var, info )
      IF (info .lt. 0) RETURN
c
c     ex-post Skewness
      CALL SKEWN(n, x, expr, var, expsk, info)
      expsk = expsk
      IF (info .lt. 0) RETURN
c    
      RETURN
      END
