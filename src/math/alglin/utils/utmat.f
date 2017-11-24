c=======================================================================
c
c     Linear Algebra - Utilities Functions                   
c
c     Note: used vectorized matrix operations for dimension < 1000
c
c-----------------------------------------------------------------------
c
c     sub-programs names convention 
c     -----------------------------
c
c     first letter : operation type
c                   A : Arrangement
c                   C : Conversion
c                   D : Difference
c                   E : Extremum
c                   I : initalization
c                   J : Inversion
c                   M : Mean 
c                   N : Norm 
c                   O : multiple Operations
c                   P : Product
c                   R : Transpose
c                   S : Sum
c                   T : Trace
c                   W : who ?
c                   X : scalar product
c                   Y : copy
c                   Z : symetrization
c  (A)
c        AVOC    : sort a vector in increase order
c        AVOD    : sort a vector in decrease order
c        AVOCI   : sort a vector in increase order (index in old vector)
c        AVODI   : sort a vector in decrease order (index in old vector)
c  (C)
c        CMCML   : converting a full square matrix (n*n) in low triangular matrix (n*n)
c        CMCMS   : converting a full square matrix (n*n) in symmetric matrix (n*(n+1)/2)
c        CMCMU   : converting a full square matrix (n*n) in up triangular matrix (n*n)
c        CMDV    : converting the diagonal of a matrix(n*n) in vector(n)
c        CMDSV   : converting the diagonal of  a symmetric matrix(n*(n+1)/2) in vector(n)
c        CMLMC   : converting a low triangular matrix (n*n) in full square matrix (n*n)
c        CMSMC   : converting a symmetric matrix (n*(n+1)/2) in full square matrix (n*n)
c        CMUMC   : converting an up triangular matrix (n*n) in full square matrix (n*n)
c        CVMD    : converting a vector (n) in diagonal matrix (n*n)
c        CVMDS   : converting a vector (n) in symmetric diagonal matrix (n*(n+1)/2)
c  (D)
c        DM      : difference of 2 matrices (vectorized full)
c        DMS     : difference of 2 matrices (vectorized symmetric)
c  (E)
c        EMBOR   : matrix (general) minimum and maximum elements
c        EMMAX   : matrix (general) maximum element
c        EMMIN   : matrix (general) minimum element
c        EMDBOR  : matrix (diagonal) minimum and maximum elements
c        EMDMAX  : matrix (diagonal) maximum element
c        EMDMIN  : matrix (diagonal) minimum element
c        EVBOR   : vector minimum and maximum elements
c        EVMAX   : vector maximum element
c        EVMAXI  : vector maximum element and index
c        EVMIN   : vector minimum element
c        EVMINI  : vector minimum element and index
c  (I)
c        IMX     : matrix (general) scalar initialization
c        IMDX    : matrix (diagonal) scalar initialization
c        IMSX    : matrix (symmetric) scalar initialization
c        IVX     : vector scalar initialization
c  (J)
c        JMC     : matrix inverse, invA(n*n) <- A(n*n)^(-1)
c        JMS     : matrix (symmetric) inverse, invA(n*n) <- A(n*n)^(-1)
c  (M)
c        MCM     : matrix (general) mean on each column
c        MCMI    : matrix (general) mean on the i-th column
c        MLM     : matrix (general) mean on each row
c        MLMI    : matrix (general) mean on the j-th row
c        MM      : matrix (general) mean 
c        MMS     : matrix (symmetric) mean
c        MV      : vector mean
c  (N)
c        NDM     : Frobenius norm of the difference between two matrices (vectorized full)
c        NDMS    : Frobenius norm of the difference between two matrices (vectorized symmetric)
c        NDVL1   : L1-norm of the difference between two vectors
c        NDVL2   : L2-norm of the difference between two vectors
c        NM      : matrix Frobenius norm (vectorized full)
c        NM1     : matrix L1-norm (vectorized full)
c        NM2     : matrix L2-norm (vectorized full)
c        NMINF   : matrix Infinite-norm (vectorized full)
c        NMS     : matrix Frobenius norm (vectorized symmetric)
c        NSV     : L2-norm of the sum of two vectors
c        NV      : vector L2-norm
c        NV1     : vector L1-norm
c        NV2     : vector L2-norm square
c        NVINF   : vector Infinite-norm
c  (S)
c        SCM     : sum of the elements of the columns of a matrix
c        SEM     : sum of the elements of a matrix
c        SEMDX   : sum of the diagonal elements of a matrix with a scalar : B = A + x*I
c                  ( A square matrix(n*n), x scalar, gives B square matrix(n*n) )
c        SEV     : sum of the elements of a vector
c        SEVP    : sum of the elements of a part of a vector
c        SLM     : sum of the elements of the rows of a matrix
c        SM      : M1 + M2     -> M (general)
c        SMbis   : sum of 2 matrices ??? 
c        SMS     : M1 + M2     -> M (vectorized symmetric)
c        SMX     : M + x       -> M
c  (T)
c        TM      : Trace of a matrix (general)
c        TMS     : Trace of a matrix (symmetric)
c  (W)
c        WVEQ    : who is equal to x in a vector ? (+/- epsilon) gives a vector of indexes
c        WVGT    : who is greater than x in a vector ? gives a vector of indexes
c        WVLT    : who is lower than x in a vector ? gives a vector of indexes
c  (X)
c        XM      : scalar product of two matrices (general)
c                  (M1 matrix n*m, M2 matrix n*m, gives M1.M2= scalar)
c        XMS     : scalar product of two matrices (symmetric)
c                  (M1 sym. matrix n*m, M2 sym.matrix n*m, gives scalar)
c  (Y)
c        YCMV    : copy a column of a vectorized matrix in a vector
c        YLMV    : copy a row of a vectorized matrix in a vector
c        YMCPI   : copy a part of a square matrix in a smaller square matrix with index
c        YMCPIR  : copy a square matrix in a part of a bigger square matrix with index
c        YMP     : copy a part of a vectorized matrix in a vectorized matrix
c        YMPMP   : copy a part of a vectorized matrix in a part of a vectorized matrix
c        YVCM    : copy a vector in a column of a vectorized matrix
c        YVLM    : copy a vector in a row of a vectorized matrix
c        YVP     : copy a part of a vector in a vector
c        YVPIR   : copy a vector in a part of a bigger vector as index
c        YVPVP   : copy a part of a vector in a part of a vector
c
c=======================================================================
c
c     subroutine AVOC
c
c     Sort a vector in increase order
c
c-----------------------------------------------------------------------
      SUBROUTINE AVOC ( n, x, y )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n  : size (n>0)                                     integer
c            x  : input vector (n)                                double
c
c     OUTPUT 
c            y  : sorted vector (n)                               double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION x(*), y(*)
c
c     local variables
      INTEGER i, j, indic, nout
      DOUBLE PRECISION a
c
c-----------------------------------------------------------------------
c
      y(1) = x(1)
c
c     if n=1 exit
      IF (n .EQ. 1) RETURN
c      
      nout = 2
      indic = 1
      DO i = 2,n
         a = x(i)
c
c         rank  
          indic = 1
          DO WHILE ( (indic .LT. nout) .AND. (a .GT. y(indic)) )
             indic = indic + 1
          ENDDO
c
c         decalage
          DO j = nout,indic+1,-1
             y(j) = y(j-1)
          ENDDO
c
          y(indic) = a
          nout = nout + 1
       ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine AVOD
c
c     Sort a vector in decrease order
c
c-----------------------------------------------------------------------
      SUBROUTINE AVOD ( nvect, vecin, vecout )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vecin  : vector to sort (nvect)                      double
c
c     OUTPUT 
c            vecout : sorted vector (nvect)                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      integer nvect
      double precision vecin(*), vecout(*)
c
      integer i, j, indic, nout
      double precision a
c
c-----------------------------------------------------------------------
c
      vecout(1) = vecin(1)
c
      if (nvect.gt.1) then
c
         nout = 2
         indic = 1
         vecout(1) = vecin(1)
         do i=2,nvect
            a = vecin(i)
c
c           recherche du rang 
c
            indic = 1
            do while ( (indic.lt.nout) .and. (a.lt.vecout(indic)) )
               indic = indic + 1
            end do
c
c           decalage
c
            do j=nout,indic+1,-1
               vecout(j) = vecout(j-1)
            end do
c
            vecout(indic) = a
            nout = nout + 1
c
         end do
c
      endif
c
      return
      end
c
c=======================================================================
c
c     subroutine AVOCI
c
c     Sort a vector in increase order with index in old vector
c
c-----------------------------------------------------------------------
      SUBROUTINE AVOCI ( nvect, vecin, vecout, ind )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vecin  : vector to sort (nvect)                      double
c
c     OUTPUT 
c            vecout : sorted vector (nvect)                       double
c            index  : indices of vector elements                 integer
c                     in the old vector       vector(nvect)
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nvect
      integer ind(*)
      double precision vecin(*), vecout(*)
c
      integer i, j, indic, nout
      double precision a
c
c-----------------------------------------------------------------------
c
      ind(1) = 1
      vecout(1) = vecin(1)
c
      if (nvect.gt.1) then
c
         nout = 2
         indic = 1
         do i=2,nvect
            a = vecin(i)
c
c           recherche du rang 
c
            indic = 1
            do while ( (indic.lt.nout) .and. (a.gt.vecout(indic)) )
               indic = indic + 1
            end do
c
c           decalage
c
            do j=nout,indic+1,-1
               vecout(j) = vecout(j-1)
               ind(j) = ind(j-1)
            end do
c
            vecout(indic) = a
            ind(indic) = i
            nout = nout + 1
c
         end do
c
      endif
c
      return
      end
c
c=======================================================================
c
c     subroutine AVODI
c
c     Sort a vector in decrease order with index in old vector
c
c-----------------------------------------------------------------------
      SUBROUTINE AVODI ( nvect, vecin, vecout, index )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vecin  : vector to sort (nvect)                      double
c
c     OUTPUT 
c            vecout : sorted vector (nvect)                       double
c            index  : indices of vector elements                 integer
c                     in the old vector       vector(nvect)
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nvect
      integer index(*)
      double precision vecin(*), vecout(*)
c
      integer i, j, indic, nout
      double precision a
c
c-----------------------------------------------------------------------
c
      index(1) = 1
      vecout(1) = vecin(1)
c
      if (nvect.gt.1) then
c
         nout = 2
         indic = 1
         do i=2,nvect
            a = vecin(i)
c
c           recherche du rang 
c
            indic = 1
            do while ( (indic.lt.nout) .and. (a.lt.vecout(indic)) )
               indic = indic + 1
            end do
c
c           decalage
c
            do j=nout,indic+1,-1
               vecout(j) = vecout(j-1)
               index(j) = index(j-1)
            end do
c
            vecout(indic) = a
            index(indic) = i
            nout = nout + 1
c
         end do
c
      endif
c
      return
      end
c
c=======================================================================
c
c     subroutine CMCML
c
c     converting a full square matrix (n*n) in low triangular matrix (n*n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMCML ( n, A, ltA )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : matrix size                                        integer
c       A   : matrix (n*n)                                        double
c
c     OUTPUT 
c       ltA : low triangular matrix (n*n)                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k
      DOUBLE PRECISION A(n,*), ltA(n,*)
c
      k = 0
      DO j = 1,n
         DO i = 1,n
            IF (i .GE. j) THEN
               ltA(i,j) = A(i,j)
            ELSE
               ltA(i,j) = 0.0
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMCMS
c
c     converting a full square matrix (n*n) in symmetric matrix (n*(n+1)/2)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMCMS ( n, A, B )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : matrix size                                          integer
c       A : matrix (n*n)                                          double
c
c     OUTPUT 
c       B : symmetric matrix (n*(n+1)/2)                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k
      DOUBLE PRECISION A(n,*), B(*)
c
      k = 0
      DO j = 1,n
         DO i = 1,j
            k = k + 1
            B(k) = A(i,j)
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMCMU
c
c     converting a full square matrix (n*n) in upper triangular matrix (n*n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMCMU ( n, A, utA )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : matrix size                                        integer
c       A   : matrix (n*n)                                        double
c
c     OUTPUT 
c       utA : upper triangular matrix (n*n)                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k
      DOUBLE PRECISION A(n,*), utA(n,*)
c
      k = 0
      DO j = 1,n
         DO i = 1,n
            IF (i .LE. j) THEN
               utA(i,j) = A(i,j)
            ELSE
               utA(i,j) = 0.
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMDV
c
c     Converting the diagonal of a matrix(n*n) in vector(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMDV ( n, A, diagA )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       A     : matrix (n*n)                                      double
c
c     OUTPUT 
c       diagA : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION A(n,*), diagA(*)
c
      DO i = 1,n
          diagA(i) = A(i,i)
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMDSV
c
c     Converting the diagonal of a symmetric matrix(n*(n+1)/2) in vector(n)
c
c
c=======================================================================
      subroutine CMDSV ( n, A, diagA )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       A     : symmetric matrix (n*(n+1)/2)                      double
c
c     OUTPUT 
c       diagA : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, k
      DOUBLE PRECISION A(*), diagA(*)
c
      k = 0
      DO i = 1,n
         k = k + i
         diagA(i) = A(k)
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMLMC
c
c     converting a low triangular matrix (n*n) in full square matrix (n*n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMLMC ( n, ltA, A )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : matrix size                                        integer
c       ltA : low triangular matrix (n*n)                         double
c
c     OUTPUT 
c       A   : matrix (n*n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j
      DOUBLE PRECISION A(n,*), ltA(n,*)
c
      DO j = 1,n
         DO i = j,n
            A(i,j) = ltA(i,j)
            IF (i .NE. j) THEN
                A(j,i) = A(i,j)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMSMC
c
c     Converting a symmetric matrix (n*(n+1)/2) in full square matrix (n*n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CMSMC ( n, symA, A)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       symA : symmetric matrix (n*(n+1)/2)                       double
c
c     OUTPUT 
c       A    : matrix (n*n)                                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k
      DOUBLE PRECISION symA(*), A(n,*)
c
      k = 0
      DO j = 1,n
         DO i = 1,j
            k = k + 1
            A(i,j) = symA(k)
            IF ( i .NE. j ) A(j,i) = symA(k)
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CMUMC
c
c     converting upper triangular matrix (n*n) to matrix (n*n)
c
c-----------------------------------------------------------------------
      subroutine CMUMC ( n, utA, A )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : matrix size                                        integer
c       utA : upper triangular matrix (n*n)                       double
c
c     OUTPUT 
c       A   : matrix (n*n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j
      DOUBLE PRECISION A(n,*), utA(n,*)
c
      DO j = 1,n
         DO i = 1,j
            A(i,j) = utA(i,j)
            IF (i .NE. j) THEN
              A(j,i) = A(i,j)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CVMD
c
c     Converting a vector (n) in full diagonal matrix (n*n)
c
c-----------------------------------------------------------------------
      SUBROUTINE CVMD ( n, V, A)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : vector size                                          integer
c       V : vector (n)                                            double
c
c     OUTPUT 
c       A : diagonal matrix (n*n)                                 double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION alpha, A(n,*), V(*)
      PARAMETER (alpha = 0.0)
c
      CALL PMX2 ( n, n, A, alpha)
      DO i = 1,n
         A(i,i) = V(i)
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CVMDS
c
c     Converting a vector (n) in symmetric diagonal matrix (n*(n+1)/2)
c
c-----------------------------------------------------------------------
      SUBROUTINE CVMDS ( n, V, A)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : vector size                                          integer
c       V : vector (n)                                            double
c
c     OUTPUT 
c       A : diagonal matrix (n*(n+1)/2)                           double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j, k
      DOUBLE PRECISION A(*), V(*)
c
      A(1) = V(1)
      IF ( n .GT. 1) THEN
         k = 1
         DO j = 2,n
            DO i = 1,j-1
               k = k + 1
               A(k) = 0.0
            ENDDO
            k = k + 1
            A(k) = V(j)
         ENDDO
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine DM
c
c     A - B -> C (difference of two general matrices)
c
c-----------------------------------------------------------------------
      SUBROUTINE DM ( n, m, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n : number of rows of the matrices             integer
c            m : number of columns of the matrices          integer
c            A : matrix (n*m)                                double
c            B : matrix (n*m)                                double
c
c     OUTPUT 
c            C : A - B, matrix (n*m)                         double
c
c-----------------------------------------------------------------------
c
      implicit none
c
c     i/o arguments
      INTEGER n, m
      DOUBLE PRECISION A(n,*), B(n,*), C(n,*)
c
      INTEGER i, j
c
c-----------------------------------------------------------------------
c
      DO j = 1,m
         DO i = 1,n
            C(i,j) = A(i,j) - B(i,j)
         ENDDO
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine DMS
c
c     A - B -> C (difference of two vectorized symmetric matrices)
c
c-----------------------------------------------------------------------
      SUBROUTINE DMS ( n, A, B, C )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n  : size (n>0)                                     integer
c            A  : symmetric matrix as vector(n*(n+1)/2)           double
c            B  : symmetric matrix as vector(n*(n+1)/2)           double
c
c     OUTPUT 
c            C  : A - B, symmetric matrix as vector(n*(n+1)/2)    double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION A(*), B(*), C(*)
c
c     local variables
      INTEGER k
c
c-----------------------------------------------------------------------
c
      k = n*(n+1)/2
      CALL DV ( k, A, B, C )
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMBOR
c
c     matrix minimum and maximum, min A(n,m) and max A(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE EMBOR ( n, m, A, Amin, Amax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       Amin : matrix minimum                                     double
c       Amax : matrix maximum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION Amin, Amax, A(*)
c
      p = n*m
      CALL EVBOR ( p, A, Amin, Amax )
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMDBOR
c
c     matrix minimum and maximum, min A(n,n) and max A(n,n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EMDBOR ( n, A, Amin, Amax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       A    : matrix (n*n)                                       double
c
c     OUTPUT 
c       Amin : matrix minimum                                     double
c       Amax : matrix maximum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Amin, Amax, x, A(n,*)
c
      Amin = A(1,1)
      Amax = A(1,1)
      IF (n .EQ. 1) RETURN         
      DO i = 2,n
        x = A(i, i)
        IF (x .LT. Amin) Amin = x
        IF (x .GT. Amax) Amax = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMDMAX
c
c     diagonal matrix maximum, max diag[A(n,n)]
c
c-----------------------------------------------------------------------
      SUBROUTINE EMDMAX ( n, A, Amax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       A    : matrix (n*n)                                       double
c
c     OUTPUT 
c       Amax : matrix maximum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Amax, x, A(n,*)
c
      Amax = A(1,1)
      IF (n .EQ. 1) RETURN
      DO i = 2,n
        x = A(i, i)
        IF (x .GT. Amax) Amax = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMDMIN
c
c     diagonal matrix minimum, min diag[A(n,n)]
c
c-----------------------------------------------------------------------
      SUBROUTINE EMDMIN ( n, A, Amin )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       A    : matrix (n*n)                                       double
c
c     OUTPUT 
c       Amin : matrix minimum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Amin, x, A(n,*)
c
      Amin = A(1,1)
      IF (n .EQ. 1) RETURN
      DO i = 2,n
        x = A(i,i)
        IF (x .LT. Amin) Amin = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMMAX
c
c     matrix maximum, max A(n,m), A vetctorized
c
c-----------------------------------------------------------------------
      subroutine EMMAX ( n, m, A, Amax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       Amax : matrix maximum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION Amax, A(*)
c
      p = n*m
      CALL EVMAX ( p, A, Amax )
      RETURN
      END
c
c=======================================================================
c
c     subroutine EMMIN
c
c      matrix minimum, min A(n,m), A vetctorized
c
c-----------------------------------------------------------------------
      SUBROUTINE EMMIN ( n, m, A, Amin )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       Amin : matrix minimum                                     double
c 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION Amin, A(*)
c
      p = n*m
      CALL EVMIN ( p, A, Amin )
      RETURN
      END
c
c=======================================================================
c
c     subroutine EVBOR
c
c     vector minimum and maximum, max V(n) and min V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EVBOR ( n, V, Vmin, Vmax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       Vmin : vector minimum                                     double
c       Vmax : vector maximum                                     double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Vmin, Vmax, x, V(*)
c
      Vmin = V(1)
      Vmax = V(1)
      IF (n .EQ. 1) RETURN
      DO i = 2,n
        x = V(i)
        IF (x .LT. Vmin) Vmin = x
        IF (x .GT. Vmax) Vmax = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EVMAX
c
c     vector maximum, max V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EVMAX ( n, V, Vmax )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       Vmax : vector maximum element                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Vmax, x, V(*)
c
      Vmax = V(1)
      IF (n .EQ. 1) RETURN
      DO i = 2,n
        x = V(i)
        IF (x .GT. Vmax) Vmax = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EVMAXI
c
c     vector maximum element and index, Vmax(index) = max V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EVMAXI ( n, V, Vmax, ind )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       Vmax : vector maximum element                             double
c       ind  : index of maximum                                  integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, ind, i
      DOUBLE PRECISION Vmax, x, V(*)
c
      Vmax = V(1)
      ind  = 1
      IF (n .EQ. 1) RETURN
      DO i = 2,n
         x = V(i)
         IF (x .GT. Vmax) THEN
            Vmax = x
            ind  = i
         ENDIF
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EVMIN
c
c     vector minimum, max V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EVMIN ( n, V, vmin )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       Vmin : vector minimum element                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION Vmin, x, V(*)
c
      Vmin = V(1)
      IF (n .EQ. 1) RETURN
      DO i = 2,n
        x = V(i)
        IF (x .LT. Vmin) Vmin = x
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine EVMINI
c
c     vector minimum element and index, Vmin(index) = min V(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE EVMINI ( n, V, Vmin, ind )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       Vmin : vector minimum element                             double
c       ind  : index of minimum                                  integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, ind, i
      DOUBLE PRECISION Vmin, x, V(*)
c
      Vmin = V(1)
      ind  = 1
      IF (n .EQ. 1) RETURN
      DO i = 2,n
         x = V(i)
         IF (x .LT. Vmin) THEN
            Vmin = x
            ind  = i
         ENDIF
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMX
c
c     matrix scalar initialization, A(i,j) = alpha
c
c-----------------------------------------------------------------------
      SUBROUTINE IMX ( n, p, A, alpha )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of rows (n>0)                             integer
c       p     : number of columns (p>0)                          integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT 
c       A     : matrix (n*p)                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, p, i, j
      DOUBLE PRECISION alpha, A(n,*)
c
      DO j = 1,p
         DO i = 1,n
             A(i,j) = alpha
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMDX
c
c     diagonal matrix scalar initialization, D(i,j) = 0.0  
c                                                   = alpha if i=j
c-----------------------------------------------------------------------
      SUBROUTINE IMDX ( n, D, alpha )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT 
c       D     : diagonal matrix (n*n)                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION alpha, ZERO, D(n,*)
      PARAMETER (ZERO = 0.0)
c
      CALL IMX (n, n, D, ZERO)
      DO i = 1,n
         D(i,i) = alpha
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMSX
c
c     symmetric matrix scalar initialization
c
c-----------------------------------------------------------------------
      SUBROUTINE IMSX ( n, A, alpha )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT 
c       A     : symmetric matrix (n*(n+1)/2)                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, p
      DOUBLE PRECISION alpha, A(*)
c
      p = n*(n + 1)/2
      CALL IVX ( p, A, alpha )
      RETURN
      END
c
c=======================================================================
c
c     subroutine IVX
c
c     vector scalar initialization, V(n) = alpha
c
c-----------------------------------------------------------------------
      SUBROUTINE IVX ( n, V, alpha )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : vector size                                      integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT 
c       V     : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i
      DOUBLE PRECISION alpha, V(*)
c
      DO i = 1,n
          V(i) = alpha
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine JMC
c
c     matrix inverse, invA(n*n) <- A(n*n)^(-1)
c
c-----------------------------------------------------------------------
      SUBROUTINE JMC ( n, A, iwork, dwork, invA, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                     integer
c       A      : matrix(n*n)                                      double
c
c     WORKSPACE 
c       iwork  : n                                               integer
c       dwork  : n                                                double
c
c     OUTPUT
c       invA   : inverse matrix (n*n)                             double
c       info   : diagnostic argument                             integer
c
c     CALL   
c        YM     : matrix copy
c        DGETRF : LU factorization of a general matrix (LAPACK)
c        DGETRI : inverse a LU factorized matrix (LAPACK)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, info
      DOUBLE PRECISION A(*), invA(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER lwork, ppi, pdw
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      lwork =  n
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      ppi = 1
c     ppi   : pointer for pivot indices of DGETRF, DGETRI (n)
c
c     total size of iwork array = (n)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdw = 1
c     pdw   : pointer for DGETRI workspace (n)
c
c     total size of dwork array : n
c
c----------------------------------------------------------------------
c
      CALL  YM ( n, n, A, invA )
      CALL DGETRF ( n, n, invA, n, iwork(ppi), info )
      CALL DGETRI ( n, invA, n, iwork(ppi), dwork(pdw), lwork, info )
      RETURN
      END
c
c=======================================================================
c
c     subroutine JMS
c
c     matrix (symmetric) inverse, invA(n*n) <- A(n*n)^(-1)
c
c-----------------------------------------------------------------------
      SUBROUTINE JMS ( n, A, dwork, invA, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                     integer
c       A      : matrix (n*n)                                     double
c
c     WORKSPACE 
c       dwork  : n*n                                              double
c
c     OUTPUT 
c       invA   : inverse matrix (n*n)                             double
c       info   : dignostic argument                              integer
c
c     CALL   
c       YM     : matrix copy 
c       DPOTRF : Choleski factorization (LAPACK)
c       DPOTRI : inverse a Choleski factorized matrix (LAPACK)
c       CMLMC  : low triangular matrix (n*n)-> square matrix (n*n)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, info
      DOUBLE PRECISION A(*), invA(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pda
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pda = 1
c     pda   : pointer temporary matrix (n*n) 
c
c     total size of dwork array : n*n
c
c-----------------------------------------------------------------------
c
      CALL YM ( n, n, A, dwork(pda) )
      CALL DPOTRF ( 'L', n, dwork(pda), n, info )
      CALL DPOTRI ( 'L', n, dwork(pda), n, info )
      CALL CMLMC ( n, dwork(pda), invA )
      RETURN
      END
c
c=======================================================================
c
c     subroutine MCM
c
c     mean of each column of matrix (general)
c
c-----------------------------------------------------------------------
      SUBROUTINE MCM ( n, m, A, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       mean : mean vector (m)                                    double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, i, j
      double precision sum, A(n,*), mean(*)
c
      DO j = 1,m
         sum = 0.0
         DO i = 1,n
            sum = sum + A(i,j)
         ENDDO
         mean(j) = sum / n
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine MCMI
c
c     mean of an index column of matrix (general)
c
c-----------------------------------------------------------------------
      SUBROUTINE MCMI ( n, m, A, ind, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c       ind  : column index                                      integer
c
c     OUTPUT 
c       mean : mean of the index column                           double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, ind, i
      DOUBLE PRECISION mean, sum, A(n,*)
c
      IF (ind .GT. m) RETURN
      sum = 0.0
      DO i = 1,n
         sum = sum + A(i,ind)
      ENDDO
      mean = sum / n
      RETURN
      END
c
c=======================================================================
c
c     subroutine MLM
c
c     mean of each row of matrix (general)
c
c-----------------------------------------------------------------------
      SUBROUTINE MLM ( n, m, A, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       mean : mean vector (n)                                    double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, i, j
      DOUBLE PRECISION sum, mean(*), A(n,*)
c
      DO j = 1,n
         sum = 0.0
         DO i = 1,m
            sum = sum + A(j,i)
         ENDDO
         mean(j) = sum / m
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine MLMI
c
c     mean of the index row of matrix (general)
c
c-----------------------------------------------------------------------
      SUBROUTINE MLMI ( n, m, A, ind, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c       ind  : row index                                         integer
c
c     OUTPUT 
c       mean : mean of the index row                              double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, ind, i
      DOUBLE PRECISION mean, sum, A(n,*)
c
      IF (ind .GT. n) RETURN
      sum = 0.0
      DO i = 1,m
         sum = sum + A(ind,i)
      ENDDO
      mean = sum / m
      RETURN
      END
c
c=======================================================================
c
c     subroutine MM
c
c     matrix mean
c
c-----------------------------------------------------------------------
      SUBROUTINE MM ( n, m, A, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       m    : number of columns                                 integer
c       A    : matrix (n*m)                                       double
c
c     OUTPUT 
c       mean : matrix mean                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m
      DOUBLE PRECISION mean, A(n,*)
c
      CALL SEM ( n, m, A, mean )
      mean = mean / (n*m)
      RETURN
      END
c
c=======================================================================
c
c     subroutine MMS
c
c     mean of a vectorized matrix
c
c-----------------------------------------------------------------------
c
c     INPUT
c            nmat   : size of the matrix                         integer
c            mat    : matrix as vector (nlmat*ncmat)              double
c
c     OUTPUT
c            moy    : mean of the matrix                          double
c
c-----------------------------------------------------------------------
c
      subroutine MMS ( nmat, mat, moy )
c
      implicit none
c
      integer nmat
      double precision mat(*)
      double precision moy
c
      integer ssm, i, idiag, icol
      double precision som
c
c-----------------------------------------------------------------------
c
      ssm = nmat*(nmat+1)/2
      som = 0.
      icol = 1
c
      do i=1,ssm
         idiag = icol*(icol+1)/2
         if ( i.eq.idiag ) then
            som  = som + mat(i)
            icol = icol + 1
         else
            som  = som + 2.*mat(i)
         end if
      end do
c
      moy = som / (nmat*nmat)
c
      return
      end
c
c=======================================================================
c
c     subroutine MV
c
c     mean of a vector
c
c-----------------------------------------------------------------------
      SUBROUTINE MV ( n, V, mean )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       V    : vector (n)                                         double
c
c     OUTPUT 
c       mean : vector mean                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION mean, V(*)
c
      CALL SEV ( n, V, mean )
      mean = mean / n
      RETURN
      END
c
c=======================================================================
c
c     subroutine NDM
c
c     Frobenius-norm of the difference between two matrices, ||A-B||_F
c          SQRT[DIAG[(A-B)'*(A-B)]]
c
c-----------------------------------------------------------------------
      SUBROUTINE NDM ( n, m, A, B, normF )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of rows                                   integer
c       m     : number of columns                                integer
c       A     : matrix (n*m)                                      double
c       B     : matrix (n*m)                                      double
c
c     OUTPUT 
c       normF : A - B Frobenius norm                              double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, i, j
      DOUBLE PRECISION normF, x, A(n,*), B(n,*)
c
      normF = 0.
      DO i = 1,m
         DO j = 1,n
            x =  A(j,i) - B(j,i)
            normF = normF + x*x
         ENDDO
      ENDDO
      normF = SQRT(normF)
      RETURN
      END
c
c=======================================================================
c
c     subroutine NDMS
c
c     Frobenius-norm of the difference between two matrices, ||A-B||_F
c     A and B symmetrics
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of the matrix                         integer
c            mat1   : matrix 1  as vector(nmat*(nmat+1)/2)        double
c            mat2   : matrix 2  as vector(nmat*(nmat+1)/2)        double
c
c     OUTPUT 
c            result : scalar = sqrt( sum( diag(M'*M) )            double
c
c-----------------------------------------------------------------------
c
      subroutine NDMS ( nmat, mat1, mat2, normF )
c
      implicit none
      integer nmat
      double precision normF, mat1(*), mat2(*)
      integer i, im, iv1, incv, k
      double precision am
c
      normF = 0.0
      im = 1
      do i=1,nmat
         incv = 2*i - 1
         do k=1,nmat
            iv1 = im
            if ( k.le.i ) then
               iv1 = iv1 + k - 1
            else
               iv1 = iv1 + incv
               incv = incv + k
            end if
            am = mat1(iv1) - mat2(iv1)
            normF  = normF + am*am
         end do
         im = im + i
      end do
      normf = SQRT(normF)
      return
      end
c
c=======================================================================
c
c     subroutine NDVL1
c
c     L1-norm of the difference of two vectors, ||x(n) - y(n)||_L1
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect1  : vector(nvect)                               double
c            vect2  : vector(nvect)                               double
c
c     OUTPUT 
c            result :                                             double
c
c-----------------------------------------------------------------------
c
      subroutine NDVL1 ( nvect, vect1, vect2, result )
c
      implicit none
c
      integer nvect
      double precision result
      double precision vect1(*), vect2(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      result = 0.
      do i=1,nvect
         result = result + abs( vect1(i) - vect2(i) )
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine NDVL2
c
c     L2-norm of the difference of two vectors, , ||x(n) - y(n)||_L2
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect1  : vector(nvect)                               double
c            vect2  : vector(nvect)                               double
c
c     OUTPUT 
c            result : scalar = sqrt(V'*V)                         double
c                     V = vect1 - vect2
c
c-----------------------------------------------------------------------
c
      subroutine NDVL2 ( nvect, vect1, vect2, result )
c
      implicit none
c
      integer nvect
      double precision result
      double precision vect1(*), vect2(*)
c
      integer i
      double precision diff
c
c-----------------------------------------------------------------------
c
      result = 0.
      do i=1,nvect
         diff = vect1(i) - vect2(i)
         result = result + diff*diff
      end do
      result = sqrt(result)
c
      return
      end
c
c=======================================================================
c
c     subroutine NM                                          
c
c     Frobenius norm of a vectorized full matrix
c     The Frobenius norm, sometimes also called the Euclidean norm 
c    (which may cause confusion with the vector -norm which also sometimes 
c     known as the Euclidean norm), is matrix norm of an matrix defined 
c     as the square root of the sum of the absolute squares of its elements
c
c                  /   __ __           \ (1/2)
c                  |   \  \         2  |
c       || X || =  |   /_ /_  X(i,j)   |     
c                  \    i  j           /
c      
c                   _____________           
c       || X || = \/ Trace (X*X')  
c
c-----------------------------------------------------------------------
      SUBROUTINE NM ( n, p, x, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       p    : number of columns                                 integer
c       x    : matrix (n*p)                                       double
c
c     OUTPUT 
c       norm : Frobenius norm                                     double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o argument
      INTEGER n, p
      DOUBLE PRECISION norm, x(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION tmp
c
      norm = 0.0
      DO i = 1,p
         DO j = 1,n
               tmp = x(j, i)
               norm = norm + tmp*tmp
         ENDDO
      ENDDO
      norm = SQRT(norm)
      RETURN
      END
c
c=======================================================================
c
c     subroutine NM1                                         
c
c     L1-norm of a vectorized full matrix
c
c     The L1-norm is the maximum absolute column sum norm defined as: 
c                           __
c                           \
c       || X || =    max    /_ |X(i,j)|     
c                 j=1,...p   i
c
c-----------------------------------------------------------------------
      subroutine NM1 ( n, p, x, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n     : number of rows (n>0)                        integer
c            p     : number of columns (p>0)                     integer
c            x     : matrix (n*p)                                 double
c
c     OUTPUT 
c            norm  : L1-norm                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o argument
      INTEGER n, p
      DOUBLE PRECISION norm, x(n,*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION sum
c
      norm = 0.0
      DO i = 1,p
         sum = 0.0
         DO j = 1,n
               sum = sum + abs(x(j, i))
         ENDDO
         IF (sum .GT. norm) norm = sum
      ENDDO
      RETURN
      END          
c
c=======================================================================
c
c     subroutine NMINF                                       version 1.0
c
c     Infinte norm of a matrix
c     The infinite norm is the maximum absolute row sum norm defined by:   
c                           __
c                           \
c       || X || =    max    /_ |X(i,j)|     
c                 i=1,...,n  j
c
c-----------------------------------------------------------------------
      SUBROUTINE NMINF ( n, p, x, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n       : nb rows of the matrix (n>0)              integer
c            p       : nb columns of the matrix (p>0)           integer
c            x       : input matrix (n*p)                        double
c
c     OUTPUT 
c            norm    : infinite norm                             double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, p
      DOUBLE PRECISION norm, x(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION sum
c
      norm = 0.0
      DO i = 1,n
         sum = 0.0
         DO j = 1,p
            sum = sum + ABS(x(i, j))
         ENDDO
         IF (sum .GT. norm) norm = sum
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine NMS
c
c     Frobenius-norm of a vectorized symmetric matrix
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of the matrix                         integer
c            mat    : matrix  as vector(nmat*(nmat+1)/2)          double
c
c     OUTPUT 
c            result : scalar = sqrt( sum( diag(M'*M) )            double
c
c-----------------------------------------------------------------------
c
      subroutine NMS ( nmat, mat, norm )
c
      implicit none
c
      integer nmat
      double precision norm, mat(*)
      integer i, im, iv1, incv, k
      double precision am
c
      norm = 0.
      im = 1
      do i=1,nmat
         incv = 2*i - 1
         do k=1,nmat
            iv1 = im
            if ( k.le.i ) then
               iv1 = iv1 + k - 1
            else
               iv1 = iv1 + incv
               incv = incv + k
            end if
            am = mat(iv1)
            norm = norm + am*am
         end do
         im = im + i
      end do
      norm = SQRT(norm)
      return
      end
c
c=======================================================================
c
c     subroutine NSV
c
c     L2-norm of the sum of 2 vectors
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect1  : vector(nvect)                               double
c            vect2  : vector(nvect)                               double
c
c     OUTPUT 
c            result : scalar = sqrt(V'*V)                         double
c                     V = vect1 + vect2
c
c-----------------------------------------------------------------------
c
      subroutine NSV ( nvect, vect1, vect2, norm )
c
      implicit none
c
      integer nvect, i
      double precision norm, sum, vect1(*), vect2(*)
c
      norm = 0.0
      do i=1,nvect
         sum = vect1(i) + vect2(i)
         norm = norm + sum*sum
      end do
      norm = SQRT(norm)
      return
      end
c
c=======================================================================
c
c     subroutine NV
c                                     _________________
c     L2-norm of a vector, ||x|| = \/ sum of x(i)*x(i)
c                                = sqrt(x'*x)
c
c-----------------------------------------------------------------------
      subroutine NV ( n, x, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n    : size (n>0)                                   integer
c            x    : input vector (n)                              double
c
c     OUTPUT 
c            norm : L2-norm                                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION norm, x(*)
c     
c     external functions
      DOUBLE PRECISION dnrm2
      EXTERNAL dnrm2   
c
c     local variables 
      INTEGER incx
      PARAMETER( incx=1 )
c
      norm = dnrm2(n, x, incx)
      RETURN
      END 
c
c=======================================================================
c
c     subroutine NV1
c
c     L1-norm of a vector, ||x|| = sum of the elements |x(i)|
c
c-----------------------------------------------------------------------
      SUBROUTINE NV1 ( nvect, vect, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c       nvect  : size of the vectors                        integer
c       vect   : vector (nvect)                              double
c
c     OUTPUT 
c       norm : L1-norm     double
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nvect, i
      double precision norm, vect(*)
c
      norm = 0.0
      do i=1,nvect
         norm = norm + DABS(vect(i))
      end do
      return
      end      
      
c
c=======================================================================
c
c     subroutine NV2
c
c     L2-norm square of a vector
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vectors                        integer
c            vect   : vector (nvect)                              double
c
c     OUTPUT 
c            result : scalar = V'*V                     double
c
c-----------------------------------------------------------------------
c
      subroutine NV2 ( nvect, vect, norm )
c
      implicit none
c
      integer nvect, i
      double precision norm, av, vect(*)
c
      norm = 0.0
      do i=1,nvect
         av = vect(i)
         norm = norm + av*av
      end do
      return
      end
c
c=======================================================================
c
c     subroutine NVINF
c
c     Infinite-norm of a vector, ||x|| = max|x(i)|
c                                         i
c
c-----------------------------------------------------------------------
      SUBROUTINE NVINF ( n, x, norm )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n    : size (n>0)                                   integer
c            x    : input vector (n)                              double
c
c     OUTPUT 
c            norm : infinite-norm                                 double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION norm, x(*)
c
c     local parameters
      INTEGER i
      DOUBLE PRECISION tmp
c
      norm = abs(x(1))
      IF (n .EQ. 1) RETURN
      DO i = 2,n
         tmp = abs(x(i))
         IF (tmp .GT. norm) norm = tmp
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine SCM
c
c     Computes the sum of columns of a matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE SCM ( nlmat, ncmat, mat, svect )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mat    : matrix as vector (nlmat*ncmat)              double
c
c     OUTPUT 
c            svect  : sum of each column vector(ncmat)            double
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nlmat, ncmat
      double precision mat(nlmat,*)
      double precision svect(*)
c
      integer i,j
      double precision sum
c
c-----------------------------------------------------------------------
c
c     y(m) = A'(n,m)*x(n) where x(i)=1.0 i=1,...,n 
c     CALL PMTV ( n, m, A, x, y )

      do j=1,ncmat
         sum = 0.
         do i=1,nlmat
            sum = sum + mat(i,j)
         end do
         svect(j) = sum
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine SEM
c
c     sum of the elements of a matrix 
c
c-----------------------------------------------------------------------
      SUBROUTINE SEM ( n, m, A, sum )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n   : number of rows                                     integer
c       m   : number of columns                                  integer
c       A   : matrix (n*m)                                        double
c
c     OUTPUT 
c       sum : sum of the elements of A(n,m)                       double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION sum, A(*)
c
      p = n*m
      CALL SEV ( p, A, sum )
      RETURN
      END
c
c=======================================================================
c
c     subroutine SEMDX
c
c     sum matrix/scalar, B(n,n) <- A(n,n) + x*Id
c
c-----------------------------------------------------------------------
      subroutine SEMDX ( n, A, alpha, B )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : matrix size                                      integer
c       A     : matrix (n*n)                                      double
c       alpha : scalar                                            double
c
c     OUTPUT 
c       B     : matrix (n*n)                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER i, n
      DOUBLE PRECISION alpha, A(n,*), B(n,*)
c
c     copy B <- A
      CALL YM ( n, n, A, B )
      DO i = 1,n
        B(i,i) = B(i,i) + alpha
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine SEV
c
c     sum of the elements of a vector
c
c-----------------------------------------------------------------------
      SUBROUTINE SEV ( n, x, sum )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n   : vector size                                        integer
c       x   : vector (n)                                          double
c
c     OUTPUT 
c       sum : sum of the elements of x(n)                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER i, n
      DOUBLE PRECISION sum, x(*)
c
      sum = 0.0
      DO i = 1,n
         sum = sum + x(i)
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine SEVP
c
c     Computing the sum of the elements of a part of a vector
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            efirst : first element of vect to sum               integer
c            nsum   : number of element of vect to sum           integer
c
c     OUTPUT 
c            sum    : sum of the elements of the vector           double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine SEVP ( nvect, vect, efirst, nsum, sum, info )
c
      implicit none
c
      integer nvect, efirst, nsum, info
      double precision vect(*)
      double precision sum
c
      integer i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      if ( efirst+nsum-1.gt.nvect ) then
         info = -2
         return
      endif
      
      sum = 0.
      do i=1,nsum
         sum = sum + vect(efirst+i-1)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine SLM
c
c     Computes the sum of rows of a matrix
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mat    : matrix as vector (nlmat*ncmat)              double
c
c     OUTPUT 
c            svect  : sum of each row vector(nlmat)               double
c
c-----------------------------------------------------------------------
c
      subroutine SLM ( nlmat, ncmat, mat, svect )
c
      implicit none
c
      integer nlmat, ncmat
      double precision mat(nlmat,*)
      double precision svect(*)
c
      integer i,j
      double precision sum
c
c-----------------------------------------------------------------------
c
      do j=1,nlmat
         sum = 0.0
         do i=1,ncmat
            sum = sum + mat(j,i)
         end do
         svect(j) = sum
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine SM
c
c     Computing the sum of 2 matrices
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of rows of the matrices             integer
c            m      : number of columns of the matrices          integer
c            mat1   : matrix(n*m)                                 double
c            mat2   : matrix(n*m)                                 double
c
c     OUTPUT 
c            mats   : mat1 + mat2   matrix(n*m)                   double
c
c-----------------------------------------------------------------------
c
      subroutine SM ( n, m, mat1, mat2, mats )
c
      implicit none
c
      integer n, m
      double precision mat1(n,*), mat2(n,*), mats(n,*)
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      do j=1,m
         do i=1,n
            mats(i,j) = mat1(i,j) + mat2(i,j)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine SMbis
c
c     Computing the sum of 2 matrices
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of rows of the matrices             integer
c            m      : number of columns of the matrices          integer
c            mat1   : matrix(n*m)                                 double
c
c     INPUT/OUTPUT 
c            mats   : mat1 + mats   matrix(n*m)                   double
c
c-----------------------------------------------------------------------
c
      subroutine SMbis ( n, m, mat1, mats )
c
      implicit none
c
      integer n, m
      double precision mat1(n,*), mats(n,*)
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      do j=1,m
         do i=1,n
            mats(i,j) = mat1(i,j) + mats(i,j)
         end do
      end do
c
      return
      end      
c
c=======================================================================
c
c     subroutine SMS
c
c     Computing the sum of 2 vectorized symmetric matrices
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of matrices                           integer
c            vect1  : symmetric matrix as vector(nmat*(nmat+1)/2) double
c            vect2  : symmetric matrix as vector(nmat*(nmat+1)/2) double
c
c     OUTPUT 
c            vects  : vect1 + vect2
c                     symmetric matrix as vector(nmat*(nmat+1)/2) double
c
c-----------------------------------------------------------------------
c
      subroutine SMS ( nmat, vect1, vect2, vects )
c
      implicit none
c
      integer nmat
      double precision vect1(*), vect2(*), vects(*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      do j=1,nmat
         do i=1,j
            k = k + 1
            vects(k) = vect1(k) + vect2(k)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine SMX
c
c     computing the sum of a matrix with a scalar
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mat    : matrix (nlmat*ncmat)                        double
c            scal   : scalar                                      double
c
c     OUTPUT 
c            mout   : matrix (nlmat*ncmat)                        double
c
c-----------------------------------------------------------------------
c
      subroutine SMX ( nlmat, ncmat, mat, scal, mout )
c
      implicit none
c
      integer nlmat, ncmat
      double precision mat(*), mout(*)
      double precision scal
c
      integer i
c
c-----------------------------------------------------------------------
c
      do i=1,nlmat*ncmat
         mout(i) = mat(i) + scal
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine TM
c
c     Computing the trace of a full square matrix
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of matrix mat                         integer
c            mat    : square matrix(nmat*nmat)                    double
c
c     OUTPUT 
c            trace  : matrix trace                                double
c
c-----------------------------------------------------------------------
c
      subroutine TM ( nmat, mat, trace )
c
      implicit none
c
      integer nmat
      double precision mat(nmat,*), trace
c
      integer i
c
c-----------------------------------------------------------------------
c
      trace = 0.
      do i=1,nmat
         trace = trace + mat(i,i)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine TMS
c
c     Computing the trace of a symmetric matrix
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : size of matrix mat                         integer
c            vect   : symmetric matrix as vector(nmat*(nmat+1)/2) double
c
c     OUTPUT 
c            trace  : matrix trace                                double
c
c-----------------------------------------------------------------------
c
      subroutine TMS ( nmat, vect, trace )
c
      implicit none
c
      integer nmat
      double precision vect(*), trace
c
      integer i, k
c
c-----------------------------------------------------------------------
c
      trace = 0.
      k = 0
      do i=1,nmat
         k = k + i
         trace = trace + vect(k)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine WVEQ
c
c     who is equal to x in a vector ? (+/- epsilon) -> index vector
c
c-----------------------------------------------------------------------
c
c     INPUT
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            val    : value                                       double
c            eps    : tolerance                                   double
c
c     OUTPUT
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements = value(+/- epsilon),
c                     vector(nind)                               integer
c
c-----------------------------------------------------------------------
c
      subroutine WVEQ ( nvect, vect, val, eps, nind, index )
c
      implicit none
c
      integer nvect, nind
      integer index(*)
      double precision val, eps
      double precision vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
c
         if ( abs(vect(i)-val).lt.eps ) then
            nind = nind + 1
            index(nind) = i
         endif
c
      end do
c     
      return
      end
c
c=======================================================================
c
c     subroutine WVGT
c
c     who is greater than x in a vector ? gives a vector of indexes
c
c-----------------------------------------------------------------------
c
c     INPUT
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            val    : value                                       double
c
c     OUTPUT
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements > value, vector(nind)    integer
c
c-----------------------------------------------------------------------
c
      subroutine WVGT ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind
      integer index(*)
      double precision val
      double precision vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
c
         if ( vect(i).gt.val ) then
            nind = nind + 1
            index(nind) = i
         endif
c
      end do
c     
      return
      end
c
c=======================================================================
c
c     subroutine WVLT
c
c     who is lower than x in a vector ? gives a vector of indexes
c
c-----------------------------------------------------------------------
c
c     INPUT
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            val    : value                                       double
c
c     OUTPUT
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements < value, vector(nind)   integer
c
c-----------------------------------------------------------------------
c
      subroutine WVLT ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind
      integer index(*)
      double precision val
      double precision vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
c
         if ( vect(i).lt.val ) then
            nind = nind + 1
            index(nind) = i
         endif
c
      end do
c     
      return
      end
c
c=======================================================================
c
c     subroutine XM
c
c     Computing scalar product of 2 matrices
c     (M1 matrix n*m, M2 matrix n*m, gives M1.M2= scalar )
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of raws of matrices                 integer
c            ncmat  : number of columns in matrices              integer
c            mat1   : symmetric matrix 1  vector(nlmat*ncmat)     double
c            mat2   : symmetric matrix 2  vector(nlmat*ncmat)     double
c     OUTPUT
c            scapro : scalar product                              double
c
c-----------------------------------------------------------------------
c
      subroutine XM ( nlmat, ncmat, mat1, mat2, scapro )
c
      implicit none
c
      integer nlmat, ncmat
      double precision mat1(nlmat,*), mat2(nlmat,*)
      double precision scapro
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      scapro = 0.0
      do j=1,ncmat
         do i=1,nlmat
               scapro = scapro + mat1(i,j)*mat2(i,j)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine XMS
c
c     Computing scalar product of 2 symmetric matrices
c     (M1 symmetric matrix n*n, M2 symmetric matrix n*n, gives scalar )
c
c-----------------------------------------------------------------------
c
c     INPUT
c            nmat   : size of the matrix mat1 and mat2           integer
c            mat1   : symmetric matrix 1  vector(nmat*(nmat+1)/2) double
c            mat2   : symmetric matrix 2  vector(nmat*(nmat+1)/2) double
c     OUTPUT
c            scapro : scalar product                              double
c
c-----------------------------------------------------------------------
c
      subroutine XMS ( nmat, mat1, mat2, scapro )
c
      implicit none
c
      integer nmat
      double precision mat1(*), mat2(*)
      double precision scapro
c
      integer i, j, ivect
c
c-----------------------------------------------------------------------
c
      scapro = 0.0
      ivect = 0
      do i=1,nmat
         do j=1,i
            ivect = ivect + 1
            if (i.eq.j) then
               scapro = scapro + mat1(ivect)*mat2(ivect)
            else
               scapro = scapro + 2.*mat1(ivect)*mat2(ivect)
            end if
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine YCMV
c
c     Copy a column of a vectorized matrix in a vector
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            matin  : input matrix as vector(nlmat*ncmat)         double
c            ncol   : number of the column to copy  (=<ncmat)   integer
c
c     OUTPUT 
c            vecout : output vector(nlmat)                        double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine YCMV ( nlmat, ncmat, matin, ncol, vecout, info )
c
      implicit none
c
      integer nlmat, ncmat, ncol, info
      double precision matin(nlmat,*), vecout(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      if( (ncol.le.0.) .and. (ncol.gt.ncmat) ) then
        info = -2
        return
      endif
c
      do i=1,nlmat
         vecout(i) = matin(i,ncol)
      end do  
c
      return
      end
c
c=======================================================================
c
c     subroutine YLMV
c
c     Copy a row of a vectorized matrix in a vector
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            matin  : input matrix as vector(nlmat*ncmat)         double
c            nrow   : number of the row to copy  (=<nlmat)       integer
c
c     OUTPUT 
c            vecout : output vector(nlmat)                        double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine YLMV ( nlmat, ncmat, matin, nrow, vecout, info )
c
      implicit none
c
      integer nlmat, ncmat, nrow, info
      double precision matin(nlmat,*), vecout(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      if( (nrow.le.0.) .and. (nrow.gt.nlmat) ) then
        info = -2
        return
      endif 
c
      do i=1,ncmat
         vecout(i) = matin(nrow,i)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine YMCPI
c
c     Copy a part of a square matrix (nin*nin)
c     in a smaller matrix (nout*nout) as index(nout<nin)
c
c------------------------------------------------------------------------
c
c     INPUT 
c            nin    : size of the input matrix matin             integer
c            matin  : input matrix (nin*nin)                      double
c            nout   : size of the output matrix matout           integer
c            ind    : index of elements to copy vector(nout)     integer

c     OUTPUT 
c            matout : output matrix (nout*nout)                   double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine YMCPI ( nin, matin, nout, ind, matout, info )
c
      implicit none
c
      integer nin, nout, info
      integer ind(*)
      double precision matin(nin,*), matout(nout,*)
c
      integer i, j, ki, kj
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      if (nout.le.nin) then
         do j=1,nout
            kj = ind(j)
            if (kj.gt.nin) then
               info = -2
               return
            end if
            do i=1,nout
               ki = ind(i)
               if (ki.gt.nin) then
                  info = -2
                  return
               end if
               matout(i,j) = matin(ki,kj)
            end do
         end do
      else
         info = -2
      end if
c
      return
      end
c
c=======================================================================
c
c     subroutine YMCPIR
c
c     Copy square matrix (nin*nin)
c     in a part of a bigger matrix (nout*nout) as index(nin<nout)
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nin    : size of the input matrix matin             integer
c            matin  : input matrix (nin*nin)                      double
c            nout   : size of the output matrix matout           integer
c            index  : index of elements to copy vector(nin)      integer
c            matout : initial output matrix (nout*nout)           double
c
c     OUTPUT 
c            matout : modified output matrix (nout*nout)          double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine YMCPIR ( nin, matin, nout, index, matout, info )
c
      implicit none
c
      integer nin, nout, info
      integer index(*)
      double precision matin(nin,*), matout(nout,*)
c
      integer i, j, ki, kj
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      if (nout.ge.nin) then
         do j=1,nin
            kj = index(j)
            if (kj.gt.nout) then
               info = -2
               return
            end if
            do i=1,nin
               ki = index(i)
               if (ki.gt.nout) then
                  info = -2
                  return
               end if
               matout(ki,kj) = matin(i,j)
            end do
         end do
      else
         info = -2
      end if
c
      return
      end
c
c=======================================================================
c
c     subroutine YMP
c
c     Copy a part of a vectorized matrix in a vectorized matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE YMP ( nrin, ncin, matin, nrout, ncout, rfirst, cfirst,
     &                 matout, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nrin   : number of rows of the matrix matin         integer
c            ncin   : number of columns of the matrix matin      integer
c            matin  : input matrix (nrin*ncin)                    double
c            nrout  : number of rows of the matrix matout        integer
c            ncout  : number of columns of the matrix matout     integer
c            rfirst : first row of matin to copy                 integer
c            cfirst : first column of matin to copy              integer
c
c     OUTPUT 
c            matout : output matrix (nrout*ncout)                 double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nrin, ncin, nrout, ncout, rfirst, cfirst, info
      DOUBLE PRECISION matin(nrin,*), matout(nrout,*)
c
c     local variables
      INTEGER i, j, ii, jj
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF ( (rfirst+nrout-1.LE.nrin).AND.(cfirst+ncout-1.LE.ncin) ) THEN
         DO j=1,ncout
            jj=j+cfirst-1
            DO i=1,nrout
               ii=i+rfirst-1
               matout(i,j) = matin(ii,jj)
            ENDDO
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YMPMP
c
c     Copy a part of a vectorized matrix in a part of a vectorized matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE YMPMP ( nbry, nbcy, nrin, ncin, matin, rfin, cfin,
     &                   nrout, ncout, matout, rfout, cfout, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nbry   : number of rows of matin to copy            integer
c            nbcy   : number of colums of matin to copy          integer
c            nrin   : number of rows of the matrix matin         integer
c            ncin   : number of columns of the matrix matin      integer
c            matin  : input matrix (nrin*ncin)                    double
c            rfin   : first row of matin to copy                 integer
c            cfin   : first column of matin to copy              integer
c            nrout  : number of rows of the matrix matout        integer
c            ncout  : number of columns of the matrix matout     integer
c            matout : initial output matrix (nrout*ncout)         double
c            rfout  : first row of matout to write               integer
c            cfout  : first column of matout to write            integer
c
c     OUTPUT 
c            matout : output matrix (nrout*ncout)                 double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nbry, nbcy, nrin, ncin, rfin, cfin, nrout, ncout, rfout, 
     &        cfout, info
      DOUBLE PRECISION matin(nrin,*), matout(nrout,*)
c
c     local variables
      INTEGER i, j, ii, jj, ki, kj
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF( (rfin+nbry-1.le.nrin).and.(cfin+nbcy-1.le.ncin).and.
     &    (rfout+nbry-1.le.nrout).and.(cfout+nbcy-1.le.ncout) ) THEN
         DO j=1,nbcy
            jj=j+cfin-1
            kj=j+cfout-1
            DO i=1,nbry
               ii=i+rfin-1
               ki=i+rfout-1
               matout(ki,kj) = matin(ii,jj)
            ENDDO
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YVCM
c
c     Copy a vector in a column of a vectorized matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE YVCM ( nlmat, ncmat, matin, vecin, ncol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            matin  : input matrix as vector(nlmat*ncmat)         double
c            vecin  : input vector(nlmat)                         double
c            ncol   : number of the column to copy  (=<ncmat)   integer
c
c     OUTPUT 
c            matin  : input matrix with column ncol modified
c                     vector(nlmat*ncmat)                         double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nlmat, ncmat, ncol, info
      DOUBLE PRECISION matin(nlmat,*), vecin(*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF( (ncol.GT.0.) .AND. (ncol.LE.ncmat) ) THEN
         DO i = 1,nlmat
            matin(i,ncol) = vecin(i)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YVLM
c
c     Copy a vector in a row of a vectorized matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE YVLM ( nlmat, ncmat, matin, vecin, nrow, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            matin  : input matrix as vector(nlmat*ncmat)         double
c            vecin  : input vector(nlmat)                         double
c            nrow   : number of the row to copy  (=<nlmat)       integer
c
c     OUTPUT 
c            matin  : input matrix with row nrow modified
c                     vector(nlmat*ncmat)                         double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nlmat, ncmat, nrow, info
      DOUBLE PRECISION matin(nlmat,*), vecin(*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
      IF( (nrow.gt.0.) .and. (nrow.le.nlmat) ) THEN
         DO i = 1,ncmat
            matin(nrow,i) = vecin(i)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YVP
c
c     Copy a part of a vector in a vector
c
c-----------------------------------------------------------------------
      SUBROUTINE YVP ( nin, vecin, efirst, nout, vecout, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nin    : size of the input vector vecin             integer
c            vecin  : input vector (nin)                          double
c            efirst : first element of vecin to copy             integer
c            nout   : size of the output vector vecout           integer
c            vecout : initial output vector (nout)                double
c
c     OUTPUT 
c            vecout : modified output vector (nout)               double
c            info   : = 0 successful exit                        integer
c                  
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nin, nout, efirst, info
      DOUBLE PRECISION vecin(*), vecout(*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
      IF ( efirst+nout-1 .LE. nin ) THEN
         DO i = 1,nout
            vecout(i) = vecin(i+efirst-1)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YVPIR
c
c     Copy a vector in a part of a bigger vector as index
c
c-----------------------------------------------------------------------
      SUBROUTINE YVPIR ( nin, vecin, nout, ind, vecout, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nin    : size of the input vector vecin             integer
c            vecin  : input vector (nin)                          double
c            nout   : size of the output vector vecout           integer
c            ind    : index of elements to copy vector(nin)      integer
c            vecout : initial output vector (nout)                double
c
c     OUTPUT 
c            vecout : modified output vector (nout)               double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nin, nout, info, ind(*)
      DOUBLE PRECISION vecin(*), vecout(*)
c
c     local variables
      INTEGER i, ki
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF (nout .GE. nin) THEN
         DO i = 1,nin
            ki = ind(i)
            IF (ki .GT. nout) THEN
               info = -1
               RETURN
            ENDIF
            vecout(ki) = vecin(i)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine YVPIR2
c
c     Copy a vector in a part of a bigger vector as index
c
c-----------------------------------------------------------------------
      SUBROUTINE YVPIR2 ( nin, vecin, nout, ind, vecout, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            nin    : size of the input vector vecin             integer
c            vecin  : input vector (nin)                          double
c            nout   : size of the output vector vecout           integer
c            ind    : index where we start to copy               integer
c            vecout : initial output vector (nout)                double
c
c     OUTPUT 
c            vecout : modified output vector (nout)               double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nin, nout, info, ind
      DOUBLE PRECISION vecin(*), vecout(*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF (nout .GT. nin+ind-1) THEN
         DO i = 1,nin
            vecout(ind+i-1) = vecin(i)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END   
c
c=======================================================================
c
c     subroutine YVPVP
c
c     Copy a part of a vector in a part of a vector
c
c-----------------------------------------------------------------------
      subroutine YVPVP ( nbry, nin, vecin, rfin,
     &                   nout, vecout, rfout, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nbry   : number of rows of vecin to copy            integer
c            nin    : size of the input vector vecin             integer
c            vecin  : input vector (nin)                          double
c            rfin   : first row of vecin to copy                 integer
c            nout   : size of the output vector vecout           integer
c            vecout : initial output vector (nout)                double
c            rfout  : first row of vecout to copy                 integer
c
c     OUTPUT 
c            vecout : modified output vector (nout)               double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nbry, nin, nout, rfin, rfout, info
      DOUBLE PRECISION vecin(*), vecout(*)
c
c     local variables
      INTEGER i, ii, ki
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
      IF (( rfin+nbry-1.LE.nin ).AND.( rfout+nbry-1.LE.nout )) THEN
         DO i = 1,nbry
             ii=i+rfin-1
             ki=i+rfout-1
            vecout(ki) = vecin(ii)
         ENDDO
      ELSE
         info = -2
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine PVXINV
c
c     computing the product of a vector with a scalar
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            scal   : scalar                                      double
c
c     OUTPUT 
c            vout   : V*X = vector(nvect)                         double
c
c-----------------------------------------------------------------------
c
      SUBROUTINE PVXINV ( n, x, y, z )
c
      IMPLICIT NONE
c
      INTEGER n
      DOUBLE PRECISION x(*), y(*), z(*)
c
c     local variables
      INTEGER i, j, index
c
c-----------------------------------------------------------------------
c
c     sort a vector in increase order
      CALL AVOC ( n, x, y )
c
      DO i=1,n
        index = 1
        DO j=1,n
           IF  (x(i) .GT. x(j)) THEN
                index = index + 1       
           ENDIF
        ENDDO
        z(i) = y(n - index + 1)    
      ENDDO
c
      RETURN
      END
