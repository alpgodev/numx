c=======================================================================
c
c     subroutine ARETURN                                     
c
c     Compound (Geometric) Average Return
c
c     Annualized return - The annualized return of a return series is 
c     defined as the constant annual return that would result in the same 
c     compound return as the series over the time period covered by the series
c
c     Let r(1), ...,r(n) be a return series of any periodicity. 
c     Then the cumulative return (compound return) of the series is 
c         __
c     r = ||  (1 + r(i)) - 1
c        i=1,n
c
c     The constant annual return rate is:
c                           __             1/d
c     annualized_return = [ || (1 + r(i)) ]     - 1
c                          i=1,n
c 
c     where d is the number of year(s) in [1,..,n]
c 
c-----------------------------------------------------------------------
      SUBROUTINE ARETURN ( n, r, ar, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of return(s) (n > 0)                integer
c            r      : returns (n)                                 double
c            d      : year(s)                                     double
c
c     OUTPUT 
c            ar     : annualized return                           double
c            info   : diagnostic argument                        integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION d, ar
      DOUBLE PRECISION r(*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      ar   = 0.0
      d    = 12
c
c     test if d > 0
      IF (d .LE. 0.) THEN
        info = -201
        RETURN
      ENDIF 
c
c     annualized return 
      ar = 1.0
      DO i = 1,n
        IF (r(i) .LT. -1.) THEN
            info = -202
            RETURN
        ENDIF 
        ar = ar*(1. + r(i))
      ENDDO
      ar = (ar**d) - 1.0
c
      RETURN
      END
