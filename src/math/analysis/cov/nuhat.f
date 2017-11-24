c=======================================================================
c
c     subroutine NUHAT                                        
c
c     Computing useful data for the shrinkage intensity
c     For further information please refer to the specification document
c     in the directory numx/src/rne/statistic/doc/ 
c
c-----------------------------------------------------------------------      
      subroutine nuhat ( n, p, x, rhomean, cov, nuh )
c-----------------------------------------------------------------------
c
c     INPUT
c            n      : max. number of values                      integer
c            p      : number of assets (p > 1)                   integer
c            x      : asset(s) return(s) values (n*p)             double
c            rhomean: mean of the returns vector p                double    
c            cov    : covariance matrix (p*p)                     double  
c  
c     OUTPUT
c            nuh    : (p*p)- matrix containing the nu matrix      double
c
c     CALL   
c
c-----------------------------------------------------------------------
c         
      implicit none 
c
c Input/Output 
c       
      integer n, p
      double precision x(*), rhomean(*), cov(*), nuh(*)
      
c Local variables      
      double precision rhotemp1, rhotemp2, temp1, temp2, sum
      integer i, j, k
      
      do i=1,p 
        rhotemp1=rhomean(i)
        do j=1,p
            sum = 0.
            rhotemp2=rhomean(j)
            do k=1,n
                temp2 = (x((i-1)*n +k)- rhotemp1)
                temp1 = temp2 * temp2
                temp2 = temp2*(x((j-1)*n+k)- rhotemp2)           
                temp2 = temp2 - cov((i-1)*p+j)     
                temp2 = (temp1 - cov((i-1)*p+i)) * temp2
                sum = sum + temp2
            enddo
            nuh(((i-1)*p+j))= sum/n
        enddo
      enddo        
      
      return 
      end
