c=======================================================================
c
c     subroutine ETAHAT
c
c     etah= sum(i,j = 1,n, (f(i,j)-cov(i,j))^2)  
c
c-----------------------------------------------------------------------      
      subroutine etahat(p, cov, f, etah)
c-----------------------------------------------------------------------
c
c     INPUT 
c
c            p      : number of assets (p > 1)                   integer
c            cov    : covariance matrix (p*p)                     double
c            f      : shrinkage target matrix (p*p)               double
c     
c     OUTPUT 
c            etah   : scalar useful for the computation of the 
c                    shrinkage intensity                          double
c
c-----------------------------------------------------------------------
c      
      implicit none 
c     i/o variables      
      integer p 
      double precision cov(*), f(*), etah
c     local variables 
      double precision temp
      integer i, j       
      
      etah = 0.
      do i=1,p 
        do j=1,i-1
            temp = cov((i-1)*p+j)-f((i-1)*p+j)
            temp = temp * temp
            etah = etah + temp
        enddo
      enddo     
      etah = 2* etah   
c     Remark: Diagonals of f and cov are equal      
      return 
      end
