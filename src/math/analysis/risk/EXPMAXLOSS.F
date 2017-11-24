c=======================================================================
c
c     subroutine EXPMAXLOSS                                  
c
c     Maximum Loss, minimum cumulated return from the beginnin 
c     time period (over t=1,...,n).
c
c     Rate of return during the last n period (sum of n previous log-returns)
c
c     R(t,n) = ln[P(t)/P(t-n)] = ln[P(t)/P(t-1)] + ... + ln[P(t-n+1)/P(t-n)]
c            = R(t) + ... + R(t-n+1)
c
c     Return of two consecutive period: 
c
c     R(t,2) = R(t) + R(t-1)
c
c     Max. Loss  := min[sum on i R(i)]
c                := min[R(1), R(1)+R(2), ..., R(1)+...+R(n)]
c
c            R : log-returns
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPMAXLOSS ( n, r, maxloss, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values ( n > 0)                  integer
c            r      : log-returns (n)                             double
c    
c     OUTPUT 
c            maxloss  : maximum loss                              double
c            info     : diagnostic argument                      integer
c    
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION maxloss, r(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION tmp, sum
c     
c     intrinsic functions
      INTRINSIC min    
c
c-----------------------------------------------------------------------
c
c     initialization
      info    = 0
      maxloss = 0.0
      sum     = 0.0
      tmp     = 0.0
c
c     min[r(1), r(1)+r(2), ..., r(1)+...+r(n)]
      DO i = 1,n 
         sum = sum + r(i)
         tmp = min(sum, tmp)
      ENDDO
c
c     max. loss
      maxloss = tmp    
c     
      RETURN
      END
