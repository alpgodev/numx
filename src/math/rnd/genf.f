c=======================================================================
c
c     GENF                                                  
c 
c     F distribution generator
c
c     Generates a random deviate from the F (variance ratio)
c     distribution with DFN degrees of freedom in the numerator
c     and DFD degrees of freedom in the denominator.
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION genf(dfn,dfd)
c-----------------------------------------------------------------------
c
c     INPUT :
c      DFN : numerator degrees of freedom (>0)                    double
c      DFD : denominator degrees of freedom (>0)                  double
c
c     Method - Directly generates ratio of chisquare variates
c
c----------------------------------------------------------------------
c      include '../stack.h'
c      
c     scalar arguments
      DOUBLE PRECISION dfd,dfn
c
c     local scalars
      DOUBLE PRECISION xden,xnum
c     
c     changed this code to call sgamma directly
c     external functions
c      DOUBLE PRECISION genchi
c      EXTERNAL genchi
      DOUBLE PRECISION sgamma
      EXTERNAL sgamma
c     
c     executable statements
 10   xnum = 2.0*sgamma(dfn/2.0)/dfn

c      GENF = ( GENCHI( DFN ) / DFN ) / ( GENCHI( DFD ) / DFD )
      xden = 2.0*sgamma(dfd/2.0)/dfd
c     changed constant so that it will not underflow at compile time
c     while not slowing generator by using double precision or logs.
c      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 20
      IF (.NOT. (xden.LE. (1.0E-37*xnum))) GO TO 20
c      call basout(io,wte,'F: generated numbers would cause overflow')
c      WRITE (*,*) ' Numerator ',xnum,' Denominator ',xden
c      next 2 lines changed to maintain truncation of large deviates.
c      WRITE (*,*) ' GENF returning 1.0E38'
c      genf = 1.0E38
c       call basout(io,wte,' GENF returning 1.0E37')
      genf = 1.0E37
      GO TO 30

   20 genf = xnum/xden
   30 RETURN
      END
