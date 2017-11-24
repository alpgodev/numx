c=======================================================================
c
c     GENUNF                                                 
c 
c     Generate Uniform random number between LOW and HIGH
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION genunf(low, high)
c-----------------------------------------------------------------------
c
c     INPUT 
c      low  : low bound (exclusive) value to be generated         double
c      high : high bound (exclusive) value to be generated        double
c
c-----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION high, low
c     
c     external functions
      DOUBLE PRECISION ranf
      EXTERNAL ranf
c     
c     executable statements
      genunf = low + (high-low)*ranf()
c
      RETURN
      END
