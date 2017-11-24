c=======================================================================
c
c     subroutine EXPUPNB                                     
c
c     Ex-post Up Number Ratio
c
c     The Up Number Ratio is a measure of the number of periods that the
c     investment was up, when the Benchmark was up, divided by the nb of
c     periods that the Benchmark was up. The larger the ratio, the better.  
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPUPNB ( n, x, y, upnb, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       upnb  : ex-post up number ratio                           double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION upnb, x(*), y(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION t, td
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      upnb = 0.0
      t    = 0.0
      td   = 0.0
c
c     gain periods     
      DO i = 1,n
        IF ((x(i) .GE. 0).AND.(y(i) .GE. 0)) THEN
            t = t + 1.0
        ENDIF  
        IF (y(i) .GE. 0) THEN
            td = td + 1.0
        ENDIF  
      ENDDO
c      
c     computing Up Number Ratio
      IF (td .GT. 0) THEN
        upnb = t / td
      ELSE
        upnb = t
      ENDIF
c
      RETURN
      END
c
c=======================================================================
