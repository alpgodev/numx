c=======================================================================
c
c     subroutine EXPVAR                                      
c
c     Ex-post Variance
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPVAR ( n, x, var, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c          n      : size of x (n > 1)                            integer
c          x      : returns (n)                                   double
c
c     OUTPUT 
c          var    : ex-post variance                              double
c          info   : diagnostic argument                          integer
c
c     CALL   
c          VARIAN : computing the variance
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION var
      DOUBLE PRECISION x(*)
c
c     external subroutines
      EXTERNAL VARIAN
c
c-----------------------------------------------------------------------
c
c     initialization
      var  = 0.0
      info = 0
c
c     computing variance
      call VARIAN ( n, x, var, info )
c
      RETURN
      END
c
c=======================================================================
