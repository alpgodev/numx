c=======================================================================
c
c     subroutine KURTOS                                      
c
c     Excess Kurtosis (un-biased) of vector x = [x(1),...,x(n)]
c
c----------------------------------------------------------------------
c
c     INPUT :
c            n     : size of x    ( n > 3)                      integer
c            x     : vector (n)                                  double
c            meanx : mean of x                                   double
c            var   : variance of x                               double
c
c     OUTPUT :
c            kurto : excess kurtosis                             double
c            info  : diagnostic argument                        integer 
c
c-----------------------------------------------------------------------
c
      SUBROUTINE KURTOS(n, x, meanx, var, kurto, info)
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info      
      DOUBLE PRECISION x(*) 
      DOUBLE PRECISION meanx, var, kurto
c
c     local variables      
      INTEGER i
      DOUBLE PRECISION EPS
      PARAMETER ( EPS = 1.E-15 )
c     
c     intrinsic functions
      INTRINSIC dfloat 
c
c---------------------------------------------------------------------
c
c     initialization
      kurto = 0.0
c
c     test var = 0 and number of points
      IF (var .LE. EPS ) THEN
         info  = -104
         RETURN
      ENDIF
      IF (n .LE. 3 ) THEN
         info  = -2
         RETURN
      ENDIF
c
c     executable statements 
      DO i = 1,n
         kurto  = (x(i)-meanx)**4 + kurto
      ENDDO
      kurto = kurto/(var * var)
      kurto = (DFLOAT(n**2 - 2*n + 3)/((n-1)*(n-2)*(n-3))) * kurto
      kurto = kurto - 3*(DFLOAT((n-1)*(2*n-3))/(n*(n-2)*(n-3))) - 3
c
      RETURN
      END
