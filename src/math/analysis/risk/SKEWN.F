c=======================================================================
c
c     subroutine SKEWN                                       
c
c     Skewness (un-biased)
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : size of x ( n > 2)                         integer
c            x      : vector(n)                                   double
c            meanx  : mean of x                                   double
c            var    : variance of x                               double
c    
c     OUTPUT :
c            skewne : skewness                                    double
c            info   : diagnistic argument                        integer 
c    
c-----------------------------------------------------------------------
c
      SUBROUTINE SKEWN(n, x, meanx, var, skewne, info)
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, i, info
      DOUBLE PRECISION x(*)
      DOUBLE PRECISION meanx, var, skewne
c      
c     local variables
      DOUBLE PRECISION myzero
      PARAMETER ( myzero = 1.e-15 )
c     
c     intrinsic functions
      INTRINSIC dfloat    
c
c----------------------------------------------------------------------
c
c     initialization
      info   = 0
      skewne = 0.0
c
c     variance test > EPS
      IF (var .lt. myzero) THEN
         info = -104
         RETURN
      ENDIF
c
c     test sample data length 
      IF ( n .lt. 3) THEN
         info = -2
         RETURN
      ENDIF
c
c     executable statements
      DO i = 1,n
         skewne = ((x(i)-meanx)**3) + skewne
      ENDDO
      skewne = (DFLOAT(n)/((n-1)*(n-2)))*(skewne/(var**(1.5)))
c
      RETURN
      END
