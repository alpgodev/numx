c=======================================================================
c
c     GENNF                                                 
c 
c     Generate random value of Noncentral F distribution
c
c     Generates a random deviate from the  noncentral F (variance ratio)
c     distribution with DFN degrees of freedom in the numerator, and DFD
c     degrees of freedom in the denominator, and noncentrality parameter
c     XNONC.
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION gennf(dfn, dfd, xnonc)
c-----------------------------------------------------------------------
c
c     INPUT : 
c      dfn   : numerator degrees of freedom (>=1.0)                 double  
c      dfd   : denominator degrees of freedom (>0)                  double
c      xnonc : noncentrality parameter (>0)                         double
c
c     Method
c
c     Directly generates ratio of noncentral numerator chisquare variate
c     to central denominator chisquare variate.
c
c-----------------------------------------------------------------------
c
c      include '../stack.h'
c      
c     scalar arguments 
      DOUBLE PRECISION dfd, dfn, xnonc
c     
c     local scalars
      DOUBLE PRECISION xden, xnum
c     
c     external functions 
c     changed the code to call GENCHI and GENNCH directly
c      DOUBLE PRECISION genchi, gennch
c      EXTERNAL genchi, gennch
      DOUBLE PRECISION sgamma, snorm
      EXTERNAL sgamma, snorm
c     
c     executable statements
c     changed the argument checker to allow DFN = 1.0
c     in the same way as GENNCH was changed.
c      GENNF = ( GENNCH( DFN, XNONC ) / DFN ) / ( GENCHI( DFD ) / DFD )
c     changed this to call SGAMMA and SNORM directly
c     xnum = gennch(dfn,xnonc)/dfn
 10   IF (dfn.GE.1.000001) GO TO 20
c     case dfn = 1.0 - here I am treating dfn as exactly 1.0
      xnum = (snorm() + sqrt(xnonc))**2
      GO TO 30
c     case dfn > 1.0
 20   xnum = (2.0*sgamma((dfn-1.0)/2.0) + (snorm()+sqrt(xnonc))**2)/dfn
c
c     xden = genchi(dfd)/dfd
 30   xden = 2.0*sgamma(dfd/2.0)/dfd
      
c     changed constant so that it will not underflow at compile time
c     while not slowing generator by using double precision or logs.
c      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 40
      IF (.NOT. (xden.LE. (1.0E-37*xnum))) GO TO 40
      call basout(io,wte,'nf: Generated numbers would cause overflow')
c     next 2 lines changed to maintain truncation of large deviates.
c      WRITE (*,*) ' GENNF returning 1.0E38'
c      gennf = 1.0E38
       call basout(io,wte,' returning 1.0E37')
      gennf = 1.0E37
      GO TO 50

   40 gennf = xnum/xden
   50 RETURN
      END
