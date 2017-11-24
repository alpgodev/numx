c=======================================================================
c
c     subroutine EXPRET                                      
c
c     Ex-post mean return
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPRET ( n, x, expr, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : size of x (n > 0)                                 integer
c       x    : return(s) (n)                                      double
c
c     OUTPUT 
c       expr : ex-post mean return                                double
c       info : diagnostic argument                               integer
c
c     CALL 
c       MV   : computing the mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION expr
      DOUBLE PRECISION x(*)
c
c     external subroutines
      EXTERNAL MV
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     computing mean 
      CALL MV ( n, x, expr )
c
      RETURN
      END
c
c=======================================================================
