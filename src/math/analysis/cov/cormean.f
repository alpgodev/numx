c=======================================================================
c
c     subroutine CORMEAN                                        
c
c     mean of the non-diagonal elements of the correlation matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE cormean (p, cov, corrmean, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            p      : number of asset(s) (p >= 1)                integer
c            cov    : covariance matrix (p*p)                     double 
c
c     OUTPUT 
c            corrmean   : mean of the correlation matrix values   double 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, info
      DOUBLE PRECISION cov(*), corrmean    
c
c     local variables
      double precision temp
      integer i, j 
c     
c     intrinsic functions
      INTRINSIC sqrt            

c   
      if (p .gt. 1) then      
        corrmean = 0.
        DO i = 1,p-1
            DO j = i+1,p
                temp = cov((i-1)*p+j)
                temp = temp/sqrt(cov((i-1)*p+i))
                temp = temp/sqrt(cov((j-1)*p+j))
                corrmean = temp + corrmean
            ENDDO
        ENDDO
          
        corrmean = corrmean
        corrmean = 2*corrmean
        corrmean = corrmean/p
        corrmean = corrmean/(p-1)      
      else
        info = -1400
      endif
c
      RETURN
      END

