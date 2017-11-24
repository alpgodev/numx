c=======================================================================
c
c     Linear Algebra - Level 2 BLAS Routines Interfaces                 
c
c-----------------------------------------------------------------------
c
c     FUNCTION LIST
c
c     (P) 
c       PEVEV   : vector product,        z(n) <- x(n).*y(n)
c       PEVJEV  : vector product,        z(n) <- x(n)./y(n)
c       PMLV    : matrix vector product, y(n) <- A(n,n)*x(n) with A low tringular
c       PMTV    : matrix vector product, y(m) <- A'(n,m)*x(n) 
c       PMV     : matrix vector product, y(n) <- A(n,m)*x(m)
c       PMVV    : matrix vector product, z(n) <- A(n,m)*x(m) + y(n) 
c       PMVX    : matrix vector product, y(n)  <- alpha*A(n,m)*x(m)
c
c       PV      : product vector/vector, A(n,m) <- x(n)*y(m)'
c       PV1     : product vector/vector, A(n,n) <- x(n)*x(n)'
c
c=======================================================================
c
c     subroutine PEVEV
c
c     vector product, z(n) <- x(n).*y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PEVEV (n, x, y, z)
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
      INTEGER n, i
      DOUBLE PRECISION x(*), y(*), z(*)
c
      DO i = 1,n
        z(i) = x(i)*y(i)
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine PEVJEV
c
c     vector product, z(n) <- x(n)./y(n)
c
c-----------------------------------------------------------------------
      SUBROUTINE PEVJEV (n, x, y, z, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : vector size                                       integer
c       x    : vector (nvect)                                     double
c       y    : vector (nvect)                                     double
c
c     OUTPUT 
c       z    : vector (n)                                         double
c       info : diagnostic argument                               integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, info, i
      DOUBLE PRECISION EPS, x(*), y(*), z(*)
      PARAMETER ( EPS = 1.0E-30 )
c
      info = 0      
      DO i = 1,n
         IF ( ABS(y(i)) .GT. EPS ) THEN
            z(i) = x(i) / y(i)
         ELSE
            info = -1
         ENDIF
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMLV
c
c     matrix vector product, y(n) <- A(n,n)*x(n) with A low tringular
c
c-----------------------------------------------------------------------
      SUBROUTINE PMLV ( n, A, x, y )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n : matrix size                                          integer
c       A : matrix (n*n)                                          double
c       x : vector (n)                                            double
c
c     OUTPUT 
c       y : vector (n)                                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, incz
      DOUBLE PRECISION A(n,*), x(*), y(*)
      PARAMETER (incz=1)
c
c     copy y <- x
      CALL YV ( n, x, y )
c
c     y <- A*y, A low triangular 
      CALL DTRMV( 'L', 'N', 'N', n, A, n, y, incz)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMTV
c
c     matrix vector product, y(m) <- A'(n,m)*x(n) 
c
c-----------------------------------------------------------------------
      SUBROUTINE PMTV ( n, m, A, x, y )
c-----------------------------------------------------------------------
c
c     INPUT
c       n  : number of rows of the matrix A                      integer
c       m  : number of columns of the matrix A                   integer
c       A  : matrix (n*m)                                         double
c       x  : vector (n)                                           double
c
c     OUTPUT
c       y  : vector(m) = A*x                                      double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m
      DOUBLE PRECISION A(n,*), x(*), y(*)
      INTEGER incx, incy
      DOUBLE PRECISION alpha, beta
      PARAMETER( incx=1, incy=1 )
      PARAMETER( alpha=1.0d0, beta=0.0d0 )
c
      CALL DGEMV ( 'T', n, m, alpha, A, n, x, incx, beta, y, incy) 
      RETURN
      END   
c
c=======================================================================
c
c     subroutine PMV
c
c     matrix vector product, y(n) <- A(n,m)*x(m)  
c
c-----------------------------------------------------------------------
      SUBROUTINE PMV ( n, m, A, x, y )
c-----------------------------------------------------------------------
c
c     INPUT
c       n  : number of rows of the matrix A                      integer
c       m  : number of columns of the matrix A                   integer
c       A  : matrix (n*m)                                         double
c       x  : vector (m)                                           double
c
c     OUTPUT
c       y  : vector(n) = A*x (n)                                  double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m
      DOUBLE PRECISION A(n,*), x(*), y(*)
      INTEGER incx, incy
      DOUBLE PRECISION alpha, beta
      PARAMETER( incx=1, incy=1 )
      PARAMETER( alpha=1.0d0, beta=0.0d0 )
c   
      CALL DGEMV ( 'N', n, m, alpha, A, n, x, incx, beta, y, incy)   
      RETURN
      END
c
c=======================================================================
c
c     subroutine PMVV
c
c     matrix vector multiply: z(n) <- A(n,m)*x(m) + y(n) 
c
c-----------------------------------------------------------------------
      SUBROUTINE PMVV ( n, m, A, x, y, z )
c-----------------------------------------------------------------------
c
c     INPUT
c       n : number of rows of the matrix A                       integer
c       m : number of columns of the matrix A                    integer
c       A : matrix (n*m)                                          double
c       x : vector (m)                                            double
c       y : vector (n)                                            double
c
c     OUTPUT
c       z : vector (n)                                            double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m
      DOUBLE PRECISION A(n,*), x(*), y(*), z(*)
      INTEGER incx, incz
      DOUBLE PRECISION alpha, beta
      PARAMETER( incx=1, incz=1 )
      PARAMETER( alpha=1.0d0, beta=1.0d0 )
c     
c     copies z <- y
      CALL YV ( n, y, z )
c		
c     z = A*x + z
      CALL DGEMV ( 'N', n, m, alpha, A, n, x, incx, beta, z, incz)      
      RETURN
      END
c=======================================================================
c
c     subroutine PMVX                   
c
c     matrix vector product, y(n) <- alpha*A(n,m)*x(m)       
c
c-----------------------------------------------------------------------
      SUBROUTINE PMVX(n, m, A, x, alpha, y, dwork, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of rows                                   integer
c       m     : number of columns                                integer
c       A     : matrix (n*m)                                      double 
c       x     : vector (m)                                        double
c       alpha : scalar                                            double
c
c     WORKSPACE
c       dwork : m 
c
c     OUTPUT
c       y     : vector (n)                                        double 
c       info  : diagnostic argument                              integer 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o                    
      INTEGER n, m, info
      DOUBLE PRECISION alpha, A(n,*), x(*), y(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     dwork(m) <- alpha*x(m)
      CALL PVX(m, x, alpha, dwork(1))
c        
c     y(n) <- A(n,m)*dwork(m)
      CALL PMV(n, m, A, dwork(1), y)
      info = 0
      RETURN
      END         
c      
c=======================================================================
c
c     subroutine PV
c
c     product vector/vector A(n,m) <- x(n)*y(m)'
c
c-----------------------------------------------------------------------
      subroutine PV ( n, m, x, y, A )
c-----------------------------------------------------------------------
c
c     INPUT
c       n  : size of vector x                                    integer
c       m  : size of vector y                                    integer
c       x  : column vector (n)                                    double
c       y  : raw vector (m)                                       double
c
c     OUTPUT
c       A  : matrix (n,m)                                         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n, m, incx, incy
      DOUBLE PRECISION alpha, ZERO, x(*), y(*), A(n,*)
      PARAMETER (alpha=1.0, incx=1, incy=1, ZERO = 0.0 )
c
c     A <- 0
      CALL PMX2 ( n, m, A, ZERO)
c
c     A <- alpha*x*y' + A
      CALL DGER(n, m, alpha, x, incx, y, incy, A, n)
      RETURN
      END
c
c=======================================================================
c
c     subroutine PV1
c
c     product vector/vector, A(n,n) <- x(n)*x(n)'
c
c-----------------------------------------------------------------------
      SUBROUTINE PV1 ( n, x, A )
c-----------------------------------------------------------------------
c
c     INPUT
c       n : vector size                                          integer
c       x : column vector (n)                                     double
c
c     OUTPUT
c       A : matrix (n,n)                                          double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION x(*), A(n,*)
c
      CALL PV ( n, n, x, x, A )
      RETURN
      END      
