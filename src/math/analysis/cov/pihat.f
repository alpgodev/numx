c=======================================================================
c
c     subroutine PIHAT                                        
c
c     Comptation of variables useful for the shrinkage intensity
c     For further informations please refer to the specification document
c     in the directory numx/src/rne/statistic/doc/ 
c
c-----------------------------------------------------------------------
      SUBROUTINE PIHAT (n, p, x, rhomean, cov, pih, rhohat) 
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values                            integer
c            p      : number of assets (p > 1)                    integer
c            x      : asset(s) return(s) values (n*p)              double
c            rhomean: mean of the returns over the period p-vector double
c            cov    : covariance matrix (p*p)                      double 
c     
c     OUTPUT 
c            pih    : see specification documentation              double
c            rhohat : see specification documentation              double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p
      DOUBLE PRECISION x(*), rhomean(*), cov(*), pih, rhohat
      
c
c     local variables
      double precision rhotemp1, rhotemp2, sum, temp,covtemp
      integer i, j, k
      
c
c-----------------------------------------------------------------------
c

      pih = 0.
      rhohat = 0. 
c
c     pih is symetric so 
c      
      do i=1,p
        rhotemp1 = rhomean(i)
        do j=1,i
            sum = 0.
            if (i .eq. j) then
                covtemp = cov((i-1)*p+i)
                 do k=1,n
                    temp = (x((i-1)*n +k)- rhotemp1)
                    temp = temp*temp
c                temp = (x(i +(k-1)*p)- rhotemp1)
c                temp = temp*(x(j+(k-1)*p)- rhotemp2)        
                    temp = temp - covtemp
                    temp = temp * temp
                    sum = sum + temp
                enddo
                sum = sum/n
                pih = pih + sum
                rhohat = rhohat + sum
            else                 
                rhotemp2 = rhomean(j)
                covtemp = cov((i-1)*p+j)
                do k=1,n
                    temp = (x((i-1)*n +k)- rhotemp1)
                    temp = temp*(x((j-1)*n+k)- rhotemp2)  
                    temp = temp - covtemp
                    temp = temp * temp
                    sum = sum + temp
                enddo
                sum = sum/n
                pih = pih + 2*sum
            endif 
        enddo
      enddo    
c
      RETURN
      END
