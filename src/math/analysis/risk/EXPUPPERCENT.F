c=======================================================================
c
c     subroutine EXPUPPERCENT                                     
c
c     Ex-post Up Percentage Ratio
c
c     The Up Percentage Ratio is a measure of the number of periods that the
c     investment outperformed the Benchmark when the Benchmark was up, 
c     divided by the nb of periods that the Benchmark was up. The larger 
c     the ratio, the better.  
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPUPPERCENT ( n, x, y, upp, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       upp   : ex-post up percentage ratio                       double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION upp, x(*), y(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION t, td
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      upp  = 0.0
      t    = 0.0
      td   = 0.0
c
c     loss periods     
      DO i = 1,n
        IF ((x(i) .GE. y(i)).AND.(y(i) .GE. 0.0)) THEN
            t = t + 1.0
        ENDIF  
        IF (y(i) .GE. 0.0) THEN
            td = td + 1.0
        ENDIF  
      ENDDO
c      
c     computing Up Percentage Ratio
      IF (td .GT. 0) THEN
        upp = t / td
      ELSE
        upp = t
      ENDIF  
c
      RETURN
      END
c
c=======================================================================
