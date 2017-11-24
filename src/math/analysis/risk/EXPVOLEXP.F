c=======================================================================
c
c     subroutine EXPVOLEXP                                   
c
c     Ex-post Exponentially Weighted Volatility
c
c----------------------------------------------------------------------
      SUBROUTINE EXPVOLEXP ( n, x, lambda, vol, info)
c----------------------------------------------------------------------
c
c     INPUT 
c        n      : size of x (n > 1)                             integer
c        x      : value(s) (n)                                   double
c        lambda : decay factor (0 <= lambda <= 1)                double  
c
c     OUTPUT 
c        vol   : ex-post exp. wieghted volatility                double
c        info  : diagnostic argument                            integer
c
c     CALL   
c        VOLAT : volatility (standard deviation)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION vol, lambda, x(*)
c
c     local parameters
      INTEGER i
      DOUBLE PRECISION tmp, mean, var      
c
c-----------------------------------------------------------------------
c
c     initialization
      vol  = 0.0
      info = 0
c
c     sum [lambda^(i-1)] 
      IF (lambda .EQ. 0.0) THEN
        tmp = 1.0
      ELSE
        tmp = 0.0
        DO i = 1,n
            tmp = tmp + (lambda**(i-1))
        ENDDO
        tmp = 1./(tmp - 1.)      
      ENDIF  
c
c     mean
      mean = 0.
      DO i = 1,n
         mean = mean + x(i)
      ENDDO
      mean = mean/n
c
c     variance
      var = 0.
      DO i = 1,n
        var = var + tmp*(lambda**(i-1))*((x(i) - mean)**2)
      ENDDO 
c
c     volatility (standard deviation)
      vol = sqrt(var)
c
      RETURN
      END
c
c=======================================================================
