c=======================================================================
c
c     Linear Algebra - Level 3 BLAS Routines Interfaces               
c
c-----------------------------------------------------------------------
c
c  (A)
c        ACMI    : matrix sorted by index-columns
c        ALMI    : matrix sorted by index-rows
c  (O)
c        OMCDMCT : C(n,n) <- M(n,n)*D(n)*M'(n,n)
c        OMCTQMC : C(n,n) <- M'(n,n)*Q(n,n)*M(n,n)
c        OMDMT   : C(n,n) <- M(n,m)*D(m)*M'(n,m)
c        OVTMCV  : alpha  <- V(n)'*M(n,n)*V(n)
c  (P)
c        PM      : matrix product, C(m,n) <- A(m,k)*B(k,n)
c        PMC     : matrix product, C(n,n) <- A(n,n)*B(n,n)
c        PMDV    : matrix product, C(n,n) <- A(n,n)*diag(D(n))
c        PMMT    : matrix product, C(m,n) <- A(m,k)*B'(n,k)
c        PMRM    : matrix product, C(n,n) <- A(n,m)*A'(n,m)
c        PMS     : matrix product, C(n,n) <- A(n*(n+1)/2)*B(n*(n+1)/2), A, B symmetric matrices
c        PMS2    : matrix product, C(n*(n+1)/2) <- A(n*(n+1)/2)*A(n*(n+1)/2), A, C symmetric
c        PMTM    : matrix product, C(n,m) <- A'(k,n)*B(k,m)
c        PRMM    : matrix product, C(m,m) <- A'(n,m)*A(n,m)
c  (R)
c        RM      : matrix transpose M'(m,n) <- M(n,m)
c
c=======================================================================
c
c     subroutine ACMI
c
c     Matrix sorted by index-columns
c
c-----------------------------------------------------------------------
      SUBROUTINE ACMI ( n, m, A, ind, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : number of rows                                     integer
c       m   : number of columns                                  integer
c       A   : input matrix (n*m)                                  double
c       ind : indices of columns to sort (m)                     integer
c
c     OUTPUT 
c       C   : sorted matrix (n*m)                                 double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, ind(*), i, j
      DOUBLE PRECISION A(n,*), C(n,*)
      DO j = 1,m
         DO i = 1,n
            C(i,j) = A(i,ind(j))
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine ALMI
c
c     Matrix sorted by index-rows 
c
c-----------------------------------------------------------------------
      SUBROUTINE ALMI ( n, m, A, ind, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : number of rows                                     integer
c       m   : number of columns                                  integer
c       A   : input matrix (n*m)                                  double
c       ind : indices of rows to sort (n)                        integer
c
c     OUTPUT 
c       C   : sorted matrix (n*m)                                 double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, ind(*), i, j
      DOUBLE PRECISION A(n,*), C(m,*)
      DO j = 1,m
         DO i = 1,n
            C(i,j) = A(ind(i),j)
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine OMCDMCT
c
c     C(n,n) <- M(n,n)*D(n)*M(n,n)', D vector n of diagonal matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE OMCDMCT ( n, M, D, mwork, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size of the matrix                              integer
c       M      : input square matrix M (n*n)                      double
c       D      : vector of the diagonal matrix D (n)              double
c
c     WORKSPACE 
c       mwork  : matrix (n*n)                                     double
c
c     OUTPUT 
c       C      : output matrix (n*n)                              double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n
      DOUBLE PRECISION D(*), M(n,*), mwork(n,*), C(n,*)
c
      INTEGER i, j
      DOUBLE PRECISION alpha, beta
      PARAMETER ( alpha = 1.D0, beta = 0.D0)
c
c-----------------------------------------------------------------------
c
c     mat(n,n) = M(n,n)*D(n) (D : diagonal matrix )
      DO j = 1,n
         DO i = 1,n
            mwork(i,j) = M(i,j)*D(j)
         ENDDO
      ENDDO
c
c     DGEMM  performs one of the matrix-matrix operations
c     out -> (M*Q)*M'
      CALL DGEMM ('N', 'T', n, n, n, alpha, mwork, n, M, n, beta, C, n)
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine OMCTQMC
c
c     C(n,n) <- M'(n,n)*Q(n,n)*M(n,n)
c
c-----------------------------------------------------------------------
      SUBROUTINE OMCTQMC ( n, M, Q, mwork, C )
c-----------------------------------------------------------------------
c     
      IMPLICIT NONE
c
      INTEGER n
      DOUBLE PRECISION Q(*), M(n,*), mwork(n,*), C(n,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER ( alpha = 1.D0, beta = 0.D0)
c     DGEMM  performs one of the matrix-matrix operations
c     out -> M'Q      
      CALL DGEMM ('T', 'N', n, n, n, alpha, M, n, Q, n, beta, mwork, n)
c     out -> (M'Q)*M   
      CALL DGEMM ('N', 'N', n, n, n, alpha, mwork, n, M, n, beta, C, n)

      RETURN
      END
c
c=======================================================================
c
c     subroutine OMDMT
c
c     C(n,n) <- A(n,m)*D(m)*A'(n,m), D vector m of diagonal matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE OMDMT ( n, m, A, D, mwork, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size 1 of the matrix A (number of rows)         integer
c       m      : size 2 of the matrix A (number of columns)      integer
c       A      : input square matrix A (n*m)                      double
c       D      : vector of the diagonal matrix D (m)              double
c
c     WORKSPACE 
c       mwork  : matrix (n*m)                                     double
c
c     OUTPUT 
c       C      : output square matrix (n*n)                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n, m
      DOUBLE PRECISION D(*), A(n,*), mwork(n,*), C(n,*)
c
      INTEGER i, j
      DOUBLE PRECISION alpha, beta
      PARAMETER ( alpha = 1.D0, beta = 0.D0)
c
c-----------------------------------------------------------------------
c
c     mat = A*D ( A=matin, D=diag : diagonal matrix ) 
c
      DO j = 1,m
         DO i = 1,n
            mwork(i,j) = A(i,j)*D(j)
         ENDDO
      ENDDO   
c      
c     DGEMM  performs one of the matrix-matrix operations
c     out -> (A*D)*A'
      CALL DGEMM ('N', 'T', n, n, m, alpha, mwork, n, A, n, beta, C, n)
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine OVTMCV
c
c     a (scalar) = V'(n)*M(n,n)*V(n) = (M'(n,n)*V(n))'*V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE OVTMCV ( n, M, V, a )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : size of the matrix and vector                      integer
c       M   : matrix as vector (n*n)                              double
c       V   : vector (n)                                          double
c
c     OUTPUT 
c       a   : V'*M*V                                              double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION a, M(n,*), V(*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION sum
c
c-----------------------------------------------------------------------
c
      a = 0.0
      DO i = 1,n
         sum = 0.0
         DO j = 1,n
            sum = sum + M(i,j)*V(j)
         ENDDO
         a = a + V(i)*sum
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine PM
c
c     General matrix product: C(m,n) <-  A(m,k)*B(k,n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PM ( m, k, n, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       m : size 1 of the matrix 1 (number of rows)              integer
c       k : size 2 of the matrix 1 (number of columns)           integer
c           = size 1 of the matrix 2
c       n : size 2 of the matrix 2 (number of columns)           integer
c       A : matrix as vector(m*k)                                 double
c       B : matrix as vector(k*n)                                 double
c
c     OUTPUT 
c       C : A*B matrix as vector(m*n)                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, k
      DOUBLE PRECISION A(m,*), B(k,*), C(m,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER (alpha = 1.D0, beta = 0.D0)
c      
c     C = A*B 
      CALL DGEMM ('N', 'N', m, n, k, alpha, A, m, B, k, beta, C, m) 
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMC
c
c     matrix product, C(n,n) <-  A(n,n)*B(n,n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMC ( n, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT
c       n : size of matrix                                       integer
c       A : matrix as vector(n*n)                                 double
c       B : matrix as vector(n*n)                                 double
c
c     OUTPUT
c       C : A*B matrix as vector(n*n)                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION A(n,*), B(n,*), C(n,*)
c
      CALL PM ( n, n, n, A, B, C )
      RETURN
      END
c=======================================================================  
c   
c     subroutine PMDV                          
c
c     matrix product, C(n,n) <- A(n,n)*diag(D(n))
c
c-----------------------------------------------------------------------
      SUBROUTINE PMDV(n, D, A, dwork, C)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       D     : vector (n)                                        double
c       A     : matrix (n*n)                                      double     
c   
c     WORKSPACE
c       dwork : n*n                                               double
c
c     OUTPUT 
c       C     : matrix (n*n)                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c      
c     arguments i/o     
      INTEGER n
      DOUBLE PRECISION D(*), A(*), C(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     transforming a vector in a diagonalized matrix
      CALL CVMD ( n, D, dwork(1))
c
c     B(n,n) <- A(n,n)*D(n,n)
      CALL PMC(n, A, dwork(1), C)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMMT
c
c     matrix product, C(n,m) = A(n,k)*B'(m,k)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMMT ( n, k, m, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : size 1 of the matrix 1 (number of rows)              integer
c       k : size 2 of the matrix 1 (number of columns)           integer
c       m : size 1 of the matrix 2 (number of rows)              integer
c       A : matrix as vector(n*k)                                 double
c       B : matrix as vector(m*k)                                 double
c
c     OUTPUT 
c       C : A*B' matrix as vector(n*m)                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, k
      DOUBLE PRECISION A(n,*), B(m,*), C(n,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER (alpha = 1.D0, beta = 0.D0)
c
c     C(m,n) = A(m,k)*B'(n,k) 
      CALL DGEMM ('N', 'T', n, m, k, alpha, A, n, B, m, beta, C, n)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMRM
c
c     matrix product, C(n,n) <- A(n,m)*A'(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMRM ( n, m, A, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : size 1 of the matrix 1 (number of rows)         integer
c       m      : size 2 of the matrix 1 (number of columns)      integer
c       A      : matrix (n*m)                                     double
c
c     OUTPUT 
c       C  : matrix (n*n)                                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n, m
      DOUBLE PRECISION A(n,*), C(n,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER ( alpha = 1.D0, beta = 0.D0)
c
c     C(n,n) = A(n,m)*A(n,m)' 
      CALL DGEMM ('N', 'T', n, n, m, alpha, A, m, A, m, beta, C, n)
      RETURN
      END            
c
c=======================================================================
c
c     subroutine PMS
c
c     matrix product, C(n,n) <- A(n*(n+1)/2)*B(n*(n+1)/2), A, B symmetric matrices
c
c-----------------------------------------------------------------------
      SUBROUTINE PMS ( n, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT
c       n : matrix size                                          integer
c       A : symmetric matrix (n*(n+1)/2)                          double
c       B : symmetric matrix (n*(n+1)/2)                          double
c
c     OUTPUT
c       C : matrix (n*n)                                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k, iv1, iv2, im, jm, incv1, incv2
      DOUBLE PRECISION sum, A(*), B(*), C(n,*)
      im = 1
      DO i = 1,n
         jm = 1
         DO j = 1,n
            sum = 0.0
            incv1 = 2*i - 1
            incv2 = 2*j - 1
            DO k = 1,n
               iv1 = im
               IF (k .LE. i) THEN
                  iv1 = iv1 + k - 1
               ELSE
                  iv1 = iv1 + incv1
                  incv1 = incv1 + k
               ENDIF
               iv2 = jm
               IF (k .LE. j) THEN
                  iv2 = iv2 + k - 1
               ELSE
                  iv2 = iv2 + incv2
                  incv2 = incv2 + k
               ENDIF
               sum  = sum + A(iv1) * B(iv2)
            ENDDO
            C(i,j) = sum
            jm = jm + j
         ENDDO
         im = im + i
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMS2
c
c     matrix product, C(n*(n+1)/2) <- A(n*(n+1)/2)*A(n*(n+1)/2) 
c     A and C symmetric
c
c-----------------------------------------------------------------------
      SUBROUTINE PMS2 ( n, A, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : matrix size                                          integer
c       A : symmetric matrix (n*(n+1)/2)                          double
c     OUTPUT 
c       C : symmetric matrix (n*(n+1)/2)                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k, k2, iv1, iv2, im, jm, incv1, incv2
      DOUBLE PRECISION sum, A(*), C(*)
      k2 = 1
      im = 1
      DO i = 1,n
         jm = 1
         DO j = 1,i
            sum = 0.0
            incv1 = 2*i - 1
            incv2 = 2*j - 1
            DO k = 1,n
               iv1 = im
               IF (k .LE. i) THEN
                  iv1 = iv1 + k - 1
               ELSE
                  iv1 = iv1 + incv1
                  incv1 = incv1 + k
               ENDIF
               iv2 = jm
               IF (k .LE. j) THEN
                  iv2 = iv2 + k - 1
               ELSE
                  iv2 = iv2 + incv2
                  incv2 = incv2 + k
               ENDIF
               sum  = sum + A(iv1) * A(iv2)
            ENDDO
            C(k2) = sum
            k2 = k2 + 1
            jm = jm + j
         ENDDO
         im = im + i
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMTM
c
c     matrix product, C(n,m) <- A'(k,n)*B(k,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMTM ( k, n, m, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c       k : size 1 of the matrix 1 (number of rows)              integer
c       n : size 2 of the matrix 1 (number of columns)           integer
c       m : size 2 of the matrix 2 (number of columns)           integer
c       A : matrix as vector(k*n)                                 double
c       B : matrix as vector(k*m)                                 double
c
c     OUTPUT 
c       C  : A'*B, matrix as vector(n*m)                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER k, n, m
      DOUBLE PRECISION A(k,*), B(k,*), C(n,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER (alpha = 1.D0, beta = 0.D0)
c
c     C(n,m) = A'(k,n)*B(k,m)
      CALL DGEMM ('T', 'N', n, m, k, alpha, A, k, B, k, beta, C, n)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PRMM
c
c     matrix product, C(m,m) <- A'(n,m)*A(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE PRMM ( n, m, A, C )
c-----------------------------------------------------------------------      
c
c     INPUT
c            n      : size 1 of the matrix 1 (number of rows)    integer
c            m      : size 2 of the matrix 1 (number of columns) integer
c                     = size 2 of the matrix 2
c            A      : matrix (n*m)                                double
c
c     OUTPUT
c            C  : A'*A matrix (m*m)                               double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m
      DOUBLE PRECISION A(n,*), C(m,*)
      DOUBLE PRECISION alpha, beta
      PARAMETER ( alpha = 1.D0, beta = 0.D0)
c
c     C(m,m) = A'(n,m)*A(n,m) 
      CALL DGEMM ('T', 'N', m, m, n, alpha, A, n, A, n, beta, C, m)
      RETURN
      END            
c      
c=======================================================================
c
c     subroutine RM
c
c     matrix transpose, A'(m*n) <- A(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE RM ( n, m, A, tA )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : number of rows                                     integer
c       m  : number of columns                                   integer
c       A  : matrix (n*m)                                         double
c
c     OUTPUT 
c       tA : transpose matrix (m*n)                               double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, i, j
      DOUBLE PRECISION A(n,*), tA(m,*)
      DO j = 1,m
         DO i = 1,n
            tA(j,i) = A(i,j)
         ENDDO
      ENDDO
      RETURN
      END
