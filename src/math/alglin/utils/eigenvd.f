c=======================================================================
c
c     subroutine EIGENVD                                    
c
c     Singular Value Decomposition   A = U * Diag(S) * VT
c
c-----------------------------------------------------------------------
      SUBROUTINE eigenvd (n, Q, dwork, D, U, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : integer, size of matrix Q. 
c       Q : (n,m)-double matrix.
c
c     WORKSPACE
c       dwork: double workspace of size n*n + 6*n + 4*n 
c
c     OUTPUT
c       D: double (n)-vector of eigenvalues.
c       U: double (n,n)-matrix U.
c       info: integer, 
c
c     CALL
c       YM, DGEEV (LAPACK)
c------------------------------------------------------------------------
c
      implicit none
c
c     i/o arguments      
      INTEGER n, info
      DOUBLE PRECISION  Q(n,*), dwork(*), U(*), D(*)
c      
c     local variables
      INTEGER lwork, ldvl, pdatmp, pdwork, pdwi, pdvl 
c
      ldvl  = 1
      lwork = 4*n + 4*n
c
c     double workspace
      pdatmp = 1
c     pdatemp needs n*n  
      pdwi = pdatmp + n*n
c     pdwi needs n  
      pdvl = pdwi + n
c     pdwl needs n  
      pdwork = pdvl + n
c     dwork needds n*n + 6*n + 4*n 
c
c     matrix copy 
      CALL YM(n, n, Q, dwork(pdatmp))
c
c     DGEEV computes for an N-by-N real nonsymmetric matrix
c     the eigenvalues and eigenvectors    
      CALL DGEEV('N', 'V', n, dwork(pdatmp), n, D, dwork(pdwi), 
     &           dwork(pdvl), ldvl, U, n, dwork(pdwork), lwork, info)
      RETURN
      END
