c=======================================================================
c
c     subroutine EXPGPERCENT                                     
c
c     Ex-post Gain Percentage Ratio
c
c     The Gain Percentage Ratio is a measure of the number of periods that the
c     investment was up divided by the number of periods that the Benchmark 
c     was up. The larger the ratio, the better.  
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPGPERCENT ( n, x, y, gp, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       gp    : ex-post gain percentage ratio                     double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION gp, x(*), y(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION t, td
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      gp   = 0.0
      t    = 0.0
      td   = 0.0
c
c     loss periods     
      DO i = 1,n
        IF (x(i) .GE. 0.0)  THEN
            t = t + 1.0
        ENDIF  
        IF (y(i) .GE. 0.0) THEN
            td = td + 1.0
        ENDIF  
      ENDDO
c      
c     computing Up Percentage Ratio
      IF (td .GT. 0) THEN
        gp = t / td
      ELSE
        gp = t
      ENDIF  
c
      RETURN
      END
