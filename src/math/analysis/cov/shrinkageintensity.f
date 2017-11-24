c=======================================================================
c
c     subroutine SHRINKAGEINTENSITY                           
c
c     Shrinkage intensity in the Ledoit-Wolf formula
c
c-----------------------------------------------------------------------
      SUBROUTINE shrinkageintensity ( n, p, X, cov, f, corrmean, dwork,
     &                                deltah )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of asset return values(p >= 1)       integer
c            p      : number of assets (p > 1)                    integer
c            X      : value(s) (n*p)                              double
c            cov    : covariance matrix (p*p)                     double
c            f      : shrinkage target matrix (p*p)               double
c            corrmean: mean of the non-diagonal elements of 
c                            the correlation matrix               double  
c     WORKSPACE 
c            dwork  : vector (p*p+p)                              double
c
c     OUTPUT  
c            deltah   : shrinkage intensity                       double  
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p
      DOUBLE PRECISION X(*), cov(*), deltah, corrmean, f(*)
c     workspace      
      double precision dwork(*)
      integer pdrhomean, pdnuh
c     local variables
      double precision rhohat, sum, temp1, temp2
      double precision etah, kappahat, pih
      double precision sqrttemp1, sqrttemp2 
      integer i, j     

c
c     external subroutines
      EXTERNAL MCM, pihat, nuhat, etahat
c     
c     intrinsic functions
      INTRINSIC sqrt            
c
c-----------------------------------------------------------------------
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdrhomean = 1
c     prho: pointer ont the mean return rho, so p
      pdnuh = pdrhomean + p
c     pdzetah : pointer for relatives rents, so ( p*p ) more
c
c     Total size of dwork array = p*p+p
c
c-----------------------------------------------------------------------
c
      CALL MCM ( n, p, X, dwork(pdrhomean))
      CALL pihat ( n, p, X, dwork(pdrhomean), cov, pih, rhohat)
      CALL nuhat ( n, p, x, dwork(pdrhomean), cov, dwork(pdnuh))
      
      sum = 0.
      do i=1,p
        do j=1,p
            if (i .ne. j) then
                temp1 = dwork(pdnuh+(i-1)*p+j-1)
                temp2 = dwork(pdnuh+(j-1)*p+i-1)
                sqrttemp1 = cov((j-1)*p+j)
                sqrttemp1 = sqrttemp1/cov((i-1)*p+i)
                sqrttemp1 = sqrt(sqrttemp1)
                sqrttemp2 = 1/sqrttemp1
                sum = sum + temp1*sqrttemp1+temp2*sqrttemp2
            endif    
        enddo
      enddo        
      sum = sum*corrmean/2
      rhohat = rhohat + sum 
      
      call etahat(p, cov, f, etah)
      
      kappahat = (pih - rhohat)/etah
      
      deltah = kappahat/n
      if (rhohat .gt. pih) then
        deltah = 0. 
      else if (deltah .gt. 1) then
        deltah = 1.
      endif
c
      RETURN
      END

