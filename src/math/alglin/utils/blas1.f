c=======================================================================
c
c     Linear Algebra - Level 1 BLAS Routines Interfaces
c
c-----------------------------------------------------------------------
c
c     FUNCTION LIST
c
c   (D)
c       DV    : vector difference, z(n) <- x(n) - y(n)
c   (P)
c       PMX   : scales a matrix by a constant, B(n,m) <- x*A(n,m)
c       PMX2  : scales a matrix by a constant, A(n,m) <- x*A(n,m)
c       PVX   : scales a vector by a constant, y(n) <- a*x(n)
c       PVX2  : scales a vector by a constant, x(n) <- a*x(n)
c   (S)
c       SV    : vector sum, z(n) <- x(n) + y(n)
c       SVVX  : vector sum, z(n) <- x(n) + alpha*y(n)
c       SVVX2 : vector sum, x(n) <- x(n) + alpha*y(n)
c       SVX   : vector sum, y(n) <- x(n) + alpha ************** (NO BLAS)
c       SVXVX : vector sum, z(n) <- alpha*[x(n) + y(n)]
c       SVXVY : vector sum, z(n) <- alpha*x(n) + beta*y(n)
c   (X)
c       XV    : dot product (scalar product), alpha <- x(n)'*y(n)
c   (Y)
c       YM    : matrix copy, B(n,m) <- A(n,m)
c       YV    : vector copy, y(n) <- x(n)
c  (Z)
c       ZMCMC : matrix symmetrization, B(n,n) <- 0.5*[A(n,n) + A(n,n)']
c       ZMCMS : matrix symmetrization, B(n*(n+1)/2) <- 0.5*[A(n,n) + A(n,n)']
c
c=======================================================================
c
c     subroutine DV
c
c     vector difference, z(n) <- x(n) - y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE DV ( n, x, y, z )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : vector size (n>0)                                    integer
c       x : vector (n)                                            double
c       y : vector (n)                                            double
c
c     OUTPUT 
c       z : vector (n)                                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION alpha, x(*), y(*), z(*)
c      
c     z(n) <- x(n) - y(n)
      alpha = -1.0
      CALL SVVX ( n, x, y, alpha, z )
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMX
c
c     scales a matrix by a constant, B(n,m) <- alpha*A(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMX ( n, m, A, alpha, B )
c-----------------------------------------------------------------------
c
c     INPUT
c       n     : number of rows of the matrix M                   integer
c       m     : number of columns of the matrix M                integer
c       A     : matrix (n,m)                                      double
c       alpha : scalar                                            double
c
c     OUTPUT
c       B     : matrix (n*m)                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION alpha, A(*), B(*)
c
      p = n*m
      CALL PVX ( p, A, alpha, B )
      RETURN
      END           
c
c=======================================================================
c
c     subroutine PMX2
c
c     scales a matrix by a constant: A(n,m) <- alpha*A(n,m)
c
c-----------------------------------------------------------------------
      SUBROUTINE PMX2 ( n, m, A, alpha)
c-----------------------------------------------------------------------      
c
c     INPUT
c       n     : number of rows                                   integer
c       m     : number of columns                                integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT
c       A     : matrix (n*m)                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, p
      DOUBLE PRECISION alpha, A(*)
c
      p = n*m
      CALL PVX2 ( p, A, alpha )
      RETURN
      END
c      
c=======================================================================
c
c     subroutine PVX
c
c      scales a vector by a constant, y(n) <- alpha*x(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PVX ( n, x, alpha, y )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : size of the vector                               integer
c       x     : vector (n)                                        double
c       alpha : scalar                                            double
c
c     OUTPUT
c       y     : vector(n)                                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx
      DOUBLE PRECISION alpha, x(*), y(*)
      PARAMETER (incx = 1)
c
c     copies a vector, v, to a vector, vout
      CALL YV ( n, x, y )
c
c     scales a vector by a constant
      CALL DSCAL( n, alpha, y, incx)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PVX2
c
c     scales a vector by a constant, x(n) <- alpha*x(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PVX2 ( n, x, alpha )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : vector size                                      integer
c       alpha : scalar                                            double
c
c     INPUT/OUTPUT 
c       x     : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx
      DOUBLE PRECISION alpha, x(*)
      PARAMETER (incx = 1)
c
c     scales a vector by a constant
      CALL DSCAL( n, alpha, x, incx)
      RETURN
      END     
c      
c=======================================================================
c
c     subroutine SV
c
c     vector sum, z(n) <- x(n) + y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE SV ( n, x, y, z )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : vector size                                          integer
c       x : vector (n)                                            double
c       y : vector (n)                                            double
c
c     OUTPUT 
c       z : vector (n)                                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx, incz
      DOUBLE PRECISION alpha, x(*), y(*), z(*)
      PARAMETER (alpha=1.0, incx=1, incz=1)
c
c     copy z <- y
      CALL YV ( n, y, z )
c
c     z <- alpha*x + z
      CALL DAXPY(n, alpha, x, incx, z, incz)
      RETURN
      END
c
c=======================================================================
c
c     subroutine SVVX
c
c     vector sum, z(n) <- x(n) + alpha*y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE SVVX ( n, x, y, alpha, z )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : vector size                                      integer
c       x     : vector (n)                                        double
c       y     : vector (n)                                        double
c       alpha : scalar                                            double
c
c     OUTPUT 
c       z     : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx, incz
      DOUBLE PRECISION alpha, x(*), y(*), z(*)
      PARAMETER (incx=1, incz=1)
c
c     copy z <- x
      CALL YV ( n, x, z )
c
c     z <- alpha*y + z
      CALL DAXPY(n, alpha, y, incx, z, incz)
      RETURN
      END      
c
c=======================================================================
c
c     subroutine SVVX2
c
c     vector sum, x(n) <- x(n) + alpha*y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE SVVX2 ( n, x, y, alpha)
c-----------------------------------------------------------------------     
c
c     INPUT 
c       n     : vector size                                      integer
c       x     : input/output vector (n)                           double
c       y     : vector (n)                                        double
c       alpha : scalar                                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx, incy
      DOUBLE PRECISION alpha, x(*), y(*)
      PARAMETER (incx=1, incy=1)
c
c     x <- alpha*y + x
      CALL DAXPY(n, alpha, y, incy, x, incx)
      RETURN
      END
c
c=======================================================================
c
c     subroutine SVX
c
c     vector sum, y(n) <- x(n) + alpha 
c
c-----------------------------------------------------------------------
      SUBROUTINE SVX ( n, x, alpha, y )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n     : vector size                                      integer
c       x     : vector (n)                                        double
c       alpha : scalar                                            double
c
c     OUTPUT 
c       y     : vector (n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER i, n
      DOUBLE PRECISION alpha, x(*), y(*)
c
      DO i = 1,n
        y(i) = x(i) + alpha
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine SVXVX
c
c     vector sum, z(n) <- alpha*[x(n) + y(n)]
c
c-----------------------------------------------------------------------
      SUBROUTINE SVXVX ( n, alpha, x, y, z )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : vector size                                     integer
c       alpha : scalar                                           double
c       x     : vector (n)                                       double
c       y     : vector (n)                                       double
c
c     OUTPUT 
c       z     : vector(n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION alpha, x(*), y(*), z(*)
c
c     z(n) <- x(n) + y(n)
      CALL SV ( n, x, y, z )
c
c     z(n) <- alpha*z(n)
      CALL PVX2 ( n, z, alpha )
      RETURN
      END
c
c=======================================================================
c
c     subroutine SVXVY
c
c     vector sum, z(n) <- alpha*x(n) + beta*y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE SVXVY ( n, x, alpha, y, beta, z )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : vector size                                     integer
c       x     : vector (n)                                       double
c       alpha : scalar                                           double
c       y     : vector (n)                                       double
c       beta  : scalar                                           double
c
c     OUTPUT 
c       z     : vector(n)                                        double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION alpha, beta, x(*), y(*), z(*)
c
c     z(n) <- alpha*x(n)
      CALL PVX2 ( n, x, alpha )
c      
c     z(n) <- z(n) + beta*y(n)
      CALL SVVX2 ( n, z, y, beta)
      RETURN
      END
c
c=======================================================================
c
c     subroutine XV
c
c     dot product a (scalar) = x(n)'*y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE XV ( n, x, y, a )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n : vector size                                          integer
c       x : raw vector (n)                                        double
c       y : column vector (n)                                     double
c
c     OUTPUT 
c       a : dot product x(n)'*y(n)                                double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION a, x(*), y(*)
c
c     local variables
      INTEGER incx, incy
      DOUBLE PRECISION DDOT
      PARAMETER(incx=1,incy=1)
c
c     external function
      EXTERNAL DDOT      
c     
c     dot product function
      a = DDOT(n, x, incx, y, incy)
      RETURN
      END
c      
c=======================================================================
c
c     subroutine YM
c
c     matrix copy, B(n,m) <- A(n,m) 
c
c-----------------------------------------------------------------------
      SUBROUTINE YM ( n, m, X, Y )
c-----------------------------------------------------------------------
c
c     INPUT
c       n  : number of rows of the matrix mat                    integer
c       m  : number of columns of the matrix mat                 integer
c       X  : input matrix as vector(n*p)                          double
c
c     OUTPUT
c       Y : output matrix as vector(n*p)                          double
c
c     CALL
c       YV : cf. blas1.f 
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER m, n, p
      DOUBLE PRECISION X(*), Y(*)
c
c-----------------------------------------------------------------------
c
      p = n*m
      CALL YV ( p, X, Y )
      RETURN
      END
c
c=======================================================================
c
c     subroutine YV
c
c     vector copy, x(n) -> y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE YV ( n, x, y )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n  : size of the vector                                  integer
c       x  : input vector (n)                                     double
c
c     OUTPUT 
c       y  : output vector (n)                                    double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incx, incy
      DOUBLE PRECISION x(*), y(*)
      PARAMETER ( incx = 1, incy = 1)
c
c     copies a vector, x, to a vector, y
      CALL DCOPY( n, x, incx, y, incy)
      RETURN
      END
c
c=======================================================================
c
c     subroutine ZMCMC
c
c     matrix symmetrization, B(n,n) <- 0.5*[A(n,n) + A(n,n)']
c
c-----------------------------------------------------------------------
      SUBROUTINE ZMCMC ( n, A, B )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : matrix size                                          integer
c       A : matrix (n*n)                                          double
c
c     OUTPUT 
c       B : matrix (n*n)                                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, i, j
      DOUBLE PRECISION A(n,*), B(n,*)
c
      DO i = 1,n
         DO j=1,n
            B(i,j) = ( A(i,j) + A(j,i) ) * 0.5
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine ZMCMS
c
c     matrix symmetrization, B(n*(n+1)/2) <- 0.5*[A(n,n) + A(n,n)']
c
c-----------------------------------------------------------------------
      SUBROUTINE ZMCMS ( n, A, B )
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
      DO i = 1,n
         DO j = 1,i
            k = k + 1
            B(k) = ( A(i,j) + A(j,i) ) * 0.5
         ENDDO
      ENDDO
      RETURN
      END
