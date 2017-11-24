c=======================================================================
c
c     subroutine EXPUPC                                     
c
c     Ex-post Up Capture
c
c     The Up Capture Ratio is a measure of the investment's compound 
c     return when the Benchmark/Index was up divided by the Benchmark's
c     compound return when the Benchmark was up. The greater the value, 
c     the better.   
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPUPC ( n, x, y, upc, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       upc   : ex-post up capture                                double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION upc, x(*), y(*)
c
c     local variables    
      INTEGER i
      DOUBLE PRECISION t, td, EPS
      PARAMETER (EPS = 1.E-30)
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      upc  = 0.0
      t    = 1.0
      td   = 1.0
c
c     gain periods     
      DO i = 1,n
        IF ((x(i) .LE. -1.).OR.(y(i) .LE. -1.)) THEN
            info = -2
            upc  = 0.
            RETURN
        ENDIF
        IF (x(i) .GE. 0) THEN
            t = t * (1.0 + x(i))
        ENDIF  
        IF (y(i) .GE. 0) THEN
            td = td * (1.0 + y(i))
        ENDIF  
      ENDDO
      t  = t - 1.
      td = td - 1.
c      
c     computing Up Capture Ratio
      IF (td .GT. EPS) THEN
        upc = t / td
      ELSE
        upc = t
      ENDIF
c
      RETURN
      END

