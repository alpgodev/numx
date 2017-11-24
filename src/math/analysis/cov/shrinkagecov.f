c=======================================================================
c
c     subroutine SHRINKAGECOV                                        
c
c     Shrinkage cov. matrix estimator with Ledoit Wolf's method
c
c-----------------------------------------------------------------------
      SUBROUTINE shrinkagecov ( n, p, X, cov, dwork, shrinkcov, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values (n > 1)                   integer
c            p      : number of asset(s) (n >= 1)                integer
c            X      : value(s) (n*p)                              double
c            cov    : covariance matrix (p*p)                     double 
c
c     WORKSPACE
c            dwork  : vector ( p*(2*p+1) )                        double
c
c     OUTPUT
c            shrinkcov   : shrinkage covariance matrix (p*p)      double  
c            info: array of dimension 1,diagnostic argument      integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION X(*), cov(*), shrinkcov(*)
c
c     workspaces     
      DOUBLE PRECISION dwork(*)
      integer pdf, pdshrink
c
c     local variables
      integer i,j 
      double precision delta, corrmean, temp
c
c     external subroutines
      EXTERNAL shrinkageintensity, shrinkagetarget
c
c-----------------------------------------------------------------------
c
      info = 0
c
c     Double workspace 

      pdf = 1 
c
c     pointer for the shrinkage target needs p*p 
      pdshrink = pdf + p*p
c     local workspace for shrinkage intensity needs p*(p+1)
c     so p*(2*p+1)
c
            
c     Total size of dwork array = p*(2*p+1) 
c
c-----------------------------------------------------------------------
c     
      call shrinkagetarget( p, cov, dwork(pdf), corrmean, info )
      if (info .lt. 0) then
        return
      endif
      
      call shrinkageintensity( n, p, X, cov, dwork(pdf), corrmean,
     &                         dwork(pdshrink), delta )
      
      do i = 1,p
        do j= 1,i-1
            temp = delta*dwork(pdf+(i-1)*p+j-1)+(1-delta)*cov((i-1)*p+j)
            shrinkcov((i-1)*p+j) = temp
            shrinkcov((j-1)*p+i) = temp
        enddo
        shrinkcov((i-1)*p+i) = cov((i-1)*p+i)
      enddo
c
      RETURN
      END

