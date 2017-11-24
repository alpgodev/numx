c=======================================================================
c
c     subroutine EXPGMEAN                                     
c
c     Ex-post Average Gain (Gain Mean)
c
c     This is a simple average (arithmetic mean) of the periods with a
c     gain. It is calculated by summing the returns for gain periods 
c     (return >= 0) and then dividing the total by the number of gain 
c     periods.     
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPGMEAN ( n, x, m, gmean, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : size of x (n > 0)                                integer
c       x     : return(s) (n)                                     double
c
c     OUTPUT 
c       m     : number of gain periods                           integer
c       gmean : ex-post gain mean                                 double
c       info  : diagnostic argument                              integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, m, info
      DOUBLE PRECISION gmean, x(*)
c
c     local variables 
      INTEGER i    
      DOUBLE PRECISION sum
c
c-----------------------------------------------------------------------
c
c     initialization
      info  = 0
      gmean = 0.0
      m     = 0
      sum   = 0.0
c
c     gain periods     
      DO i = 1,n
        IF (x(i) .GE. 0) THEN
            m   = m + 1
            sum = sum + x(i)
        ENDIF  
      ENDDO
c      
c     computing mean
      IF (m .GT. 0) gmean = sum / m
c
      RETURN
      END

