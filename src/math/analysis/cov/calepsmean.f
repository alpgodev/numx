c=======================================================================
c
c     subroutine calepsmean                                    
c
c     Calibrates Kato sensitivity parameter with mean method
c
c-----------------------------------------------------------------------
      SUBROUTINE calepsmean(n, cov, epskgr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       p : integer, number of variables (p>=1). 
c       cov : double array of dimension p by p, covariance matrix. 
c
c     OUTPUT: 
c       epskgr : double, Kato eps-group parameter.
c       info   :   = 0 successful exit                        integer
c                        -1 ncup != ncdown   
c 
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE 
c      
c     I/O
      INTEGER n, info
      DOUBLE PRECISION cov(*), epskgr
c
c     local variable
      INTEGER i
c
c     initialisation
      info = 0
      epskgr = 0.0
c
      DO i=1,n
      	epskgr = epskgr + cov((i-1)*n+i)
      ENDDO
      epskgr = epskgr /n
c
      RETURN
      END
