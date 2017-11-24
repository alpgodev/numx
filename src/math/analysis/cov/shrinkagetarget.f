c=======================================================================
c
c     subroutine SHRINKAGETARGET                         
c
c     Shrinkage target in the Ledoit-Wolf formula
c
c-----------------------------------------------------------------------
      SUBROUTINE shrinkagetarget ( p, cov, f, corrmean, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c
c            p      : number of assets (p > 1)                    integer
c            cov    : covariance matrix (p*p)                     double
c
c     OUTPUT 
c
c            f   : shrinkage target  (p*p)          double
c            corrmean: mean of the correlation matrix values      double
c            info: array of dimension 1, diagnostic argument     integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, info
      DOUBLE PRECISION cov(*), f(*), corrmean
c
c     local variables
      integer i, j
      double precision temp
c
c     external subroutines
      external cormean
c     
c     intrinsic functions
      intrinsic sqrt
c
c-----------------------------------------------------------------------
c    
c     mean of the non-diagonal elements of the correlation matrix
      call cormean ( p, cov, corrmean, info )
      if (info .lt. 0) then
        return
      endif
c
c     Shrinkage target matrix
      DO i = 1,p
         DO j = 1,i-1
            temp = corrmean*sqrt(cov((i-1)*p+i)*cov((j-1)*p+j))
            f((i-1)*p+j) = temp
            f((j-1)*p+i) = temp
         ENDDO
         f((i-1)*p+i) = cov((i-1)*p+i)
      ENDDO
c
      RETURN
      END

