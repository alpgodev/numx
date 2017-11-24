c=======================================================================
c
c     subroutine EXPDNB                                     
c
c     Ex-post Down Number Ratio
c
c     The Down Number Ratio is a measure of the number of periods that the
c     investment was down, when the Benchmark was down, divided by the nb of
c     periods that the Benchmark was down. The smaller the ratio, the better.  
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPDNB ( n, x, y, dnb, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of period(s) (n > 0)                      integer
c       x     : portfolio's return(s) (n)                         double
c       y     : benchmark's return(s) (n)                         double
c
c     OUTPUT 
c       dnb   : ex-post down number ratio                         double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION dnb, x(*), y(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION t, td
c      
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      dnb  = 0.0
      t    = 0.0
      td   = 0.0
c
c     loss periods     
      DO i = 1,n
        IF ((x(i) .LE. -1.).OR.(y(i) .LE. -1.)) THEN
            info = -2
            dnb  = 0.
            RETURN
        ENDIF
        IF ((x(i) .LT. 0.0).AND.(y(i) .LT. 0.0)) THEN
            t = t + 1.0
        ENDIF  
        IF (y(i) .LT. 0.0) THEN
            td = td + 1.0
        ENDIF  
      ENDDO
c      
c     computing Down Number Ratio
      IF (td .GT. 0) THEN   
        dnb = t / td
      ELSE
        dnb = t
      ENDIF  
c
      RETURN
      END

