c=======================================================================
c
c     subroutine EXPDPERCENT                                     
c
c     Ex-post Down Percentage Ratio
c
c     The Down Percentage Ratio is a measure of the number of periods that the
c     investment outperformed the Benchmark when the Benchmark was down, 
c     divided by the nb of periods that the Benchmark was down. The larger 
c     the ratio, the better.  
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPDPERCENT ( n, x, y, dp, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       dp    : ex-post down percentage ratio                     double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION dp, x(*), y(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION t, td
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      dp  = 0.0
      t    = 0.0
      td   = 0.0
c
c     loss periods     
      DO i = 1,n
        IF ((x(i) .GE. y(i)).AND.(y(i) .LT. 0.0)) THEN
            t = t + 1.0
        ENDIF  
        IF (y(i) .LT. 0.0) THEN
            td = td + 1.0
        ENDIF  
      ENDDO
c      
c     computing Up Percentage Ratio
      IF (td .GT. 0) THEN
        dp = t / td
      ELSE
        dp = t
      ENDIF  
c
      RETURN
      END

