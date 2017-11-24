c=======================================================================
c
c     subroutine calepsvar                                   
c
c     This function computes a level of calibration (espilon) defined by 
c     a percentage of explained variance.
c
c-----------------------------------------------------------------------
      SUBROUTINE calepsvar ( n, cov, exvar, iwork, dwork,
     &                       epsilon, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p     : number of variables (p>=1)                       integer
c       cov   : covariance matrix (p*p)                           double
c       exvar : percentage of explained variance, in ]0,1]        double      
c
c     WORKSPACE 
c       iwork : n                                                integer 
c       dwork : 2*n*n + 17*n                                      double  
c
c     OUTPUT
c       epsilon : sensitivity parameter                           double          
c       info    : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE 
c      
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION cov(n, *), exvar, epsilon
c
c     workspace 
      DOUBLE PRECISION dwork(*), iwork(*)
c      
c     local variables
      INTEGER i, pisdp, pdsdp, pdsvd, pds
      DOUBLE PRECISION sum, psum, a
c
c-----------------------------------------------------------------------
c
c     initialisation
      info = 0
      epsilon = 0.0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      pisdp = 1
c     pisdp   : pointer for TESTSDP who needs (n)
c
c     Total size of iwork array ( n )
c
c     pointers for double precision work space  : dwork
c     ----------------------------------------------------------
      pdsdp = 1
c     pdsdp   : pointer for TESTSDP so (n*(2*n + 7) 
      pdsvd = pdsdp + ( n*( 2*n + 7 ) ) 
c     pdsvd  : pointer for SVD so (9*n)
      pds   = pdsvd + ( 9*n )
c     pds    : pointer for eigenvalues so (n)
c
c     Total size of work array = ( n*(2*n + 17) )
c
c-----------------------------------------------------------------------
c
c     test if n >= 2
       IF (n .LT. 2) THEN
        info = -2 
        RETURN
       ENDIF
c       
c     covariance matrix SDP test
      CALL TESTSDP (n, cov, iwork(pisdp), dwork(pdsdp), info)
      IF (info .NE. 0) THEN
         info = -108
         RETURN
      ENDIF       
c
c     singular value decomposition (svd) of cov. matrix
c     A = U * Diag(S) * VT
      CALL svd ( n, n, cov, dwork(pdsvd), 
     &           dwork(pdsdp), dwork(pdsdp + (n*n)), dwork(pds), info)
c      CALL eigenvd (n, cov, dwork(pdsvd), dwork(pds), dwork(pdsdp),info)
c
c     sum of eigenvalues
      sum =0.0
      DO i = 1,n
        sum = sum + dwork(pds + i - 1)
      ENDDO
      a    = exvar*sum
      psum = dwork(pds)
c
c     test if lambda(max) > exvar*sum
c     the % variance is too small (risk overestimated)
      IF (psum .GT. a) THEN
        epsilon = psum
        info    =  1301
        RETURN
      ENDIF
c
c     compute optimal espsilon
      DO i = 2,n
        psum = psum + dwork(pds + i - 1)
        IF (psum .GT. a) THEN
            epsilon = dwork(pds + i - 1)
c           test % variance is too large              
            IF (i .EQ. n) THEN 
                info = 1302
            ENDIF
            RETURN
        ENDIF
      ENDDO      
c
c     test % variance is too large
      IF (psum .LT. a) THEN
        epsilon = dwork(pds + n - 1)
        info = 1302
      ENDIF
c
      RETURN
      END
