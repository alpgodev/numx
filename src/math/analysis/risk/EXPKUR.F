c=======================================================================
c
c     subroutine EXPKUR                                      
c
c     Ex-post Kurtosis 
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPKUR ( n, x, expku, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c       n       : size of x                                      integer
c       x       : sample data (n)                                 double
c       
c     OUTPUT :
c       expku   : ex-post kurtosis                                double
c       info    : diagnostic argument                            integer                
c
c     CALL :
c       EXPVAR  : ex-post robust variance
c       EXPRET  : ex-post return
c       KURTOS  : kurtosis  
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION expku
      DOUBLE PRECISION x(*)
c
c     local variables
      DOUBLE PRECISION mexpr, var      
c
c     external subroutines
      EXTERNAL EXPRET, EXPVAR, KURTOS
c
c-----------------------------------------------------------------------
c
c     initialization
      expku = 0.0
      info  = 0
c
c     ex-post return
      CALL EXPRET ( n, x, mexpr, info )
      IF (info .lt. 0) RETURN
c
c     ex-post variance
      CALL EXPVAR ( n, x, var, info)
      IF (info .lt. 0) RETURN
c
c     ex-post kurtosis
      call KURTOS (n, x, mexpr, var, expku, info)
      IF (info .lt. 0) RETURN
c
      RETURN
      END
