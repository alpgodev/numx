c=======================================================================
c
c     subroutine calepsmedian
c
c     Calibrates Kato sensitivity parameter with median method
c
c-----------------------------------------------------------------------      
      SUBROUTINE calepsmedian (n, cov, dwork, epskgr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c        p : integer, number of variables (p>=1). 
c        cov : double array of dimension p by p, covariance matrix.
c
c     WORKSPACE
c       dwork: double workspace of size  2*n*n + 9*n  
c
c     OUTPUT
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
      DOUBLE PRECISION cov(*), epskgr, dwork(*)
      CHARACTER JOBVL, JOBVR
c
c     local variables
      INTEGER ldvl, lwork, pdcov, pdeigval, pdwi, pdvl, pdeigvect, 
     & pdwork, pdeigvalorder   
c
c     external function
      EXTERNAL YM, AVOC, MEDIAN 
c
c     initialisation
      info = 0       
        
      pdcov = 1
c     needs n*n so n*n
      pdeigval = pdcov + n*n 
c     needs n no n*n + n
      pdwi = pdeigval + n   
c     needs n no n*n + 2*n
      pdvl =   pdwi +n  
c     needs n so n*n +3*n
      pdeigvect = pdvl +n  
c     needs n*n so 2*n*n +3*n
      pdeigvect = pdvl +n   
c     needs n so 2*n*n +4*n
      pdwork = pdeigvect + n * n
c     needs 4*n so 2*n*n + 8*n
      pdeigvalorder = pdwork + 4*n
c     needs n, so 2*n*n + 9*n  
c  
      CALL YM(n, n, cov, dwork(pdcov))
c      
      JOBVL = 'N'
      JOBVR = 'N'     
      ldvl = 1
      lwork = 4*n
c     
      CALL dgeev(jobvl, jobvr, n, dwork(pdcov), n, dwork(pdeigval), 
     & dwork(pdwi), dwork(pdvl), ldvl, dwork(pdeigvect), n, 
     & dwork(pdwork), lwork, info)
      CALL AVOC(n, dwork(pdeigval), dwork(pdeigvalorder))
      CALL MEDIAN(n, dwork(pdeigvalorder), epskgr, info)
      RETURN
      END
