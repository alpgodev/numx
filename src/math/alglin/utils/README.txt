c
c     Name convention (Linear Algebra)
c
c     First letters : Kind of operations
c                   A : arrangements
c                   C : conversion
c                   D : difference
c                   E : extremums
c                   I : initalisation
c                   J : inversion
c                   M : moyenne
c                   N : norme
c                   O : operations multiples
c                   P : produit
c                   R : transpose
c                   S : somme
c                   T : trace
c                   W : who ?
c                   X : produit scalaire
c                   Y : copie
c                   Z : symetrisation
c
c     Following letters : operandes
c                   BOR: bornes
c                   CM : colonne d'une matrice
c                   D  : vecteur associe a une matrice diagonale
c                   DI : eigenvectors and eigenvalues.
c                   EM : elements d'une matrice
c                   EMD: elements diagonaux d'une matrice
c                   EQ : equal
c                   EV : elements d'un vecteur
c                   GT : greater
c                   I  : index
c                   J  : inverse
c                   LM : ligne d'une matrice
c                   LT : lower
c                   M  : matrice (n*m)
c                   MAX: maximum
c                   MC : matrice carree (n*n)
c                   MCP: matrice carree partielle (n*n)
c                   MCT: matrice carree transposee (n*n)
c                   MD : matrice diagonale (n*n)
c                   MDS: matrice diagonale symetrique (n*(n+1)/2)
c                   MIN: minimum
c                   ML : matrice triangulaire inferieure (low) (n*n)
c                   MP : matrice partielle (n*m)
c                   MS : matrice symetrique (n*(n+1)/2)
c                   MT : matrice transposee (m*n)
c                   MU : matrice triangulaire superieure (up) (n*n)
c                   OC : ordre croissant
c                   OD : ordre decroissant
c                   R  : fonction inverse (du meme nom sans"R")
c                   V  : vecteur (n)
c                   VT : vecteur transpose (n)
c                   X  : scalaire
c
c   (A)
c        ACMI    : sort the columns of a vectorized matrix as index
c        ALMI    : sort the rows of a vectorized matrix as index
c        AVOC    : sort a vector in increase order
c        AVOD    : sort a vector in decrease order
c        AVOCI   : sort a vector in increase order
c                  with index in old vector
c        AVODI   : sort a vector in decrease order
c                  with index in old vector
c   (C)
c        CMCML   : converting a full square matrix (n*n)
c                  in low triangular matrix (n*n)
c        CMCMS   : converting a full square matrix (n*n)
c                  in symmetric matrix (n*(n+1)/2)
c        CMCMU   : converting a full square matrix (n*n)
c                  in up triangular matrix (n*n)
c        CMDV    : converting the diagonal of a matrix(n*n) in vector(n)
c        CMDSV   : converting the diagonal of
c                  a symmetric matrix(n*(n+1)/2) in vector(n)
c        CMLMC   : converting a low triangular matrix (n*n)
c                  in full square matrix (n*n)
c        CMSMC   : converting a symmetric matrix (n*(n+1)/2)
c                  in full square matrix (n*n)
c        CMUMC   : converting an up triangular matrix (n*n)
c                  in full square matrix (n*n)
c        CVMD    : converting a vector (n) in diagonal matrix (n*n)
c        CVMDS   : converting a vector (n) in symmetric
c                  diagonal matrix (n*(n+1)/2)
c   (D)	
c        DM      : computing the difference of 2 matrices
c        DMS     : computing the difference of 2 vectorized symmetric
c                  matrices
c        DV      : BLAS 1 - vector difference, z(n) <- x(n) - y(n)
c   (E)
c        EMBOR   : computing the minimum and maximum of a matrix
c        EMDBOR  : computing the minimum and maximum
c                  of the diagonal of a matrix
c        EMDMAX  : computing the maximum of the diagonal of a matrix
c        EMDMIN  : computing the minimum of the diagonal of a matrix
c        EMMAX   : computing the maximum of a matrix
c        EMMIN   : computing the minimum of a matrix
c        EVBOR   : computing the minimum and maximum of a vector
c        EVMAX   : computing the maximum of a vector
c        EVMAXI  : computing the maximum of a vector and index
c        EVMIN   : computing the minimum of a vector
c        EVMINI  : computing the minimum of a vector and index
c   (I)
c        IMX     : initialization at a scalar of a matrix
c        IMDX    : initialization at a scalar of a diagonal matrix
c        IMSX    : initialization at a scalar of a vectorized symmetric
c                  matrix
c        IVX     : initialization at a scalar of a vector
c   (J)
c        JMC     : Inversion of a square matrix
c                  ( A matrix n*n, gives B = 1/A  matrix n*n )
c        JMS     : Inversion of a symmetric matrix
c                  ( A matrix n*n, gives B = 1/A  matrix n*n )
c   (M)
c        MCM     : computing the mean of each column of a full matrix
c        MCMI    : computing the mean of a column of a full matrix
c        MLM     : computing the mean of each row of a full matrix
c        MLMI    : computing the mean of a row of a full matrix
c        MM      : computing the mean of a vectorized matrix
c        MMS     : computing the mean of a vectorized symmetric matrix
c        MV      : computing the mean of a vector
c   (N) 
c        NDM     : computing the Frobenius norm of the difference
c                  between two vectorized full matrix
c        NDMS    : computing the Frobenius norm of the difference
c                  between two vectorized symmetric matrix
c        NDVL1   : Computing the L_1 norm of the difference of 2 vectors
c        NDVL2   : computing the L_2 norm of the difference of 2 vectors
c        NM      : computing the Frobenius norm of
c                  a vectorized full matrix
c        NMINF   : computing the inf norm of a matrix
c        NMS     : computing the Frobenius norm of
c                  a vectorized symmetric matrix
c        NSV     : computing the L_2 norm of the sum of 2 vectors
c        NV      : computing the L_2 norm of a vector
c        NV2     : computing the L_2 norm square of a vector
c        NVINF   : computing the inf norm of a vector
c   (O)
c        OMCDMCT : BLAS 3 - C(n,n) <- M(n,n)*D(n)*M'(n,n)
c        OMCTQMC : BLAS 3 - C(n,n) <- M'(n,n)*Q(n,n)*M(n,n)
c        OMDMT   : BLAS 3 - C(n,n) <- M(n,m)*D(m)*M'(n,m)
c        OVTMCV  : BLAS 3 - x <- V(n)'*M(n,n)*V(n)
c   (P)
c        PEVEV   : BLAS 2 - vector product, z(n) <- x(n).*y(n)
c        PEVJEV  : BLAS 2 - vector product, z(n) <- x(n)./y(n)
c        PM      : BLAS 3 - matrix product, C(m,n) <-  A(m,k)*B(k,n)
c        PMC     : BLAS 3 - matrix product, C(n,n) <-  A(n,n)*B(n,n)
c        PMLV    : BLAS 2 - matrix vector product, y(n) <- A(n,n)*x(n) with A low tringular
c        PMMT    : BLAS 3 - matrix product, C(m,n) <- A(m,k)*B'(n,k)
c        PMRM    : BLAS 3 - matrix product, C(n,n) <- A(n,m)*A'(n,m)
c        PMS     : BLAS 3 - matrix product, C(n,n) <- A(n*(n+1)/2)*B(n*(n+1)/2), A, B symmetrics
c        PMS2    : BLAS 3 - matrix product, C(n,n) <- A(n*(n+1)/2)*A(n*(n+1)/2), A symmetric
c        PMTM    : BLAS 3 - matrix product, C(n,m) <- A'(k,n)*B(k,m)
c        PMTV    : BLAS 2 - matrix vector product, y(m) <- A'(n,m)*x(n)
c        PMV     : BLAS 2 - matrix vector product, y(n) <- A(n,m)*x(m)
c        PMVV    : BLAS 2 - matrix vector product: z(n) <- A(n,m)*x(m) + y(n) 
c        PMX     : BLAS 1 - scales a matrix by a constant, B(n,m) <- x*A(n,m)
c        PMX2    : BLAS 1 - scales a matrix by a constant, A(n,m) <- x*A(n,m)
c        PRMM    : BLAS 3 - matrix product, C(m,m) <- A'(n,m)*A(n,m)
c        PV      : BLAS 2 - product vector/vector, X(n,m) <- U(n)*V(m)'
c        PV1     : BLAS 2 - product vector/vector, X(n,n) <- U(n)*U(n)'
c        PVX     : BLAS 1 - scales a vector by a constant, y(n) <- a*x(n)
c        PVX2    : BLAS 1 - scales a vector by a constant, x(n) <- a*x(n)
c   (R)
c        RM      : BLAS 3 - matrix transpose, A'(m*n) <- A(n,m)
c   (S)
c        SCM     : computing the sum of the elements of the columns of a matrix
c        SEM     : computing the sum of the elements of a matrix
c        SEMDX   : computing the sum of the diagonal elements of a
c                  matrix with a scalar : B = A + x*I
c                  ( A square matrix(n*n), x scalar,
c                    gives B square matrix(n*n) )
c        SEV     : computing the sum of the elements of a vector
c        SEVP    : computing the sum of the elements of a part of a vector
c        SLM     : computing the sum of the elements of the rows of a matrix
c        SM      : computing the sum of 2 matrices
c        SMS     : computing the sum of 2 vectorized symmetric matrices
c        SMX     : computing the sum of a matrix with a scalar
c        SV      : BLAS 1 - vector sum, z(n) <- x(n) + y(n)
c        SVVX    : BLAS 1 - vector sum, z(n) <- x(n) + alpha*y(n)
c        SVVX2   : BLAS 1 - vector sum, x(n) <- x(n) + alpha*y(n)
c        SVX     : BLAS 1 - vector sum, y(n) <- x(n) + alpha
c        SVXVX   : BLAS 1 - vector sum, z(n) <- alpha*[x(n) + y(n)]
c        SVXVY   : BLAS 1 - vector sum, z(n) <- alpha*x(n) + beta*y(n)
c   (T)
c        TM      : computing the trace of a full matrix
c        TMS     : computing the trace of a symmetric matrix
c   (W)
c        WVEQ    : who is equal to x in a vector ? (+/- epsilon)
c                  gives a vector of indexes
c        WVGT    : who is greater than x in a vector ?
c                  gives a vector of indexes
c        WVLT    : who is lower than x in a vector ?
c                  gives a vector of indexes
c   (X)
c        XM      : computing scalar product of two matrices
c                  (M1 matrix n*m, M2 matrix n*m, gives M1.M2= scalar)
c        XMS     : computing scalar product of two symmetric matrices
c                  (M1 sym. matrix n*m, M2 sym.matrix n*m, gives scalar)
c        XV      : BLAS 1 - dot product (scalar product), alpha <- x(n)'*y(n)
c   (Y)
c        YCMV    : copy a column of a vectorized matrix in a vector
c        YLMV    : copy a row of a vectorized matrix in a vector
c        YM      : BLAS 1 - matrix copy,  B(n,m) <- A(n,m) 
c        YMCPI   : copy a part of a square matrix
c                  in a smaller square matrix with index
c        YMCPIR  : copy a square matrix
c                  in a part of a bigger square matrix with index
c        YMP     : copy a part of a vectorized matrix
c                  in a vectorized matrix
c        YMPMP   : copy a part of a vectorized matrix
c                  in a part of a vectorized matrix
c        YV      : BLAS 1 - vector copy, y(n) <- x(n)
c        YVCM    : copy a vector in a column of a vectorized matrix
c        YVLM    : copy a vector in a row of a vectorized matrix
c        YVP     : copy a part of a vector in a vector
c        YVPIR   : copy a vector in a part of a bigger vector as index
c        YVPVP   : copy a part of a vector in a part of a vector
c   (Z)
c        ZMCMC   : symmetrises a vectorized full square matrix
c                  and gives a vectorized full square matrix A=(A+A')/2
c        ZMCMS   : symmetrises a vectorized full matrix
c                  and gives a vectorized symmetric matrix A=(A+A')/2
c
c=======================================================================
