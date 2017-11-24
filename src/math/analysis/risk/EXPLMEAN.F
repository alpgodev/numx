c=======================================================================
c
c     subroutine EXPLMEAN                                     
c
c     Ex-post Average Loss (Loss Mean)
c
c     This is a simple average (arithmetic mean) of the periods with a
c     loss. It is calculated by summing the returns for loss periods 
c     (return < 0) and then dividing the total by the number of loss 
c     periods.     
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPLMEAN ( n, x, m, lmean, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : size of x (n > 0)                                integer
c       x     : return(s) (n)                                     double
c
c     OUTPUT 
c       m     : number of loss period(s)                         integer
c       lmean : ex-post loss mean                                 double
c       info  : diagnostic argument                              integer
c
c     CALL 
c       EXPGMEAN : gain mean
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, m, info
      DOUBLE PRECISION lmean, x(*)
c
c     local variables   
      INTEGER i  
      DOUBLE PRECISION sum
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      lmean = 0.0
      m     = 0
      sum   = 0.0
c
c     loss periods     
      DO i = 1,n
        IF (x(i) .LT. 0.0) THEN
            m   = m + 1
            sum = sum + x(i)
        ENDIF  
      ENDDO
c      
c     computing loss mean
      IF (m .GT. 0) lmean = sum / m
c
      RETURN
      END

