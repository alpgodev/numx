c=======================================================================
c
c     subroutine SVD                                     
c
c     Singular value decomposition   A(n,m) =U(n,n)*S*V(m,m)'
c
c-----------------------------------------------------------------------
      SUBROUTINE svd ( m, n, A, dwork, U, VT, S, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       m     : number of rows                                   integer
c       n     : number of columns                                integer
c       A     : matrix (n*m)                                      double
c            
c     WORKSPACE
c       dwork : 8*min(m,n)+max(m,n)                              integer
c
c     OUTPUT
c       U     : matrix (n*n)                                      double
c       VT    : matrix (m*m)                                      double
c       S     : increased sorted singular values min(n,m)         double
c       info  : diagnostic argument                              integer
c
c------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER m, n, info
      DOUBLE PRECISION A(n,*), S(*), U(n,*), VT(m,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c      
c     local variables
      INTEGER lda, lwork, ldu, ldvt, maxmn
c      
      maxmn = MAX(n,m)      
      ldu   = maxmn
      ldvt  = maxmn
      lda   = maxmn    
      lwork = 9*maxmn  
c
c     DGESVD computes the singular value decomposition (SVD) 
c     of a real N-by-M matrix        
      CALL DGESVD( 'A', 'A', n, m, A, lda, S, U, ldu, VT, ldvt,
     $             dwork(1), lwork, info )     
      RETURN
      END
