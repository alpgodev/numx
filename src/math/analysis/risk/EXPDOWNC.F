c=======================================================================
c
c     subroutine EXPDOWNC                                     
c
c     Ex-post Down Capture
c
c     The Down Capture Ratio is a measure of the investment's compound 
c     return when the Benchmark/Index was down divided by the Benchmark's
c     compound return when the Benchmark was down. The smaller the value, 
c     the better.   
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPDOWNC ( n, x, y, downc, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of period(s) (n > 0)                     integer
c       x      : portfolio's return(s) (n)                        double
c       y      : benchmark's return(s) (n)                        double
c
c     OUTPUT 
c       downnc : ex-post down capture                             double
c       info   : diagnostic argument                             integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION downc, x(*), y(*)
c
c     local variables   
      INTEGER i
      DOUBLE PRECISION t, td, EPS
      PARAMETER (eps = 1.E-30)
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      downc = 0.0
      t     = 1.0
      td    = 1.0
c
c     loss periods     
      DO i = 1,n
        IF ((x(i) .LE. -1.).OR.(y(i) .LE. -1.)) THEN
            info  = -109
            downc = 0.
            RETURN
        ENDIF
        IF (x(i) .LT. 0.0) THEN
            t = t * (1.0 + x(i))
        ENDIF  
        IF (y(i) .LT. 0.0) THEN
            td = td * (1.0 + y(i))
        ENDIF  
      ENDDO
      t  = t - 1.
      td = td - 1.
c      
c     computing Down Capture Ratio
      IF ((td .GT. EPS).OR.(td .LT. -EPS)) THEN 
        downc = t / td
      ELSE
        downc = ABS(t)
      ENDIF  
c
      RETURN
      END
