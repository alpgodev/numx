c=======================================================================
c
c     subroutine EXPVOL                                      
c
c     Ex-post Volatilty of vector x = [x(1),...,x(t)].
c
c----------------------------------------------------------------------
      SUBROUTINE EXPVOL ( n, x, vol, info)
c----------------------------------------------------------------------
c
c     INPUT 
c        n     : size of x (n > 1)                              integer
c        x     : value(s) (n)                                    double
c
c     OUTPUT 
c        vol   : ex-post volatility                             double
c        info  : diagnostic argument                            integer
c
c     CALL   
c        VOLAT : volatility (standard deviation)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION vol
      DOUBLE PRECISION x(*)
c
c     external subroutines
      EXTERNAL VOLAT
c
c-----------------------------------------------------------------------
c
c     initialization
      vol  = 0.0
      info = 0
c
c     computing volatility
      CALL VOLAT ( n, x, vol, info )
c
      RETURN
      END
c
c=======================================================================
