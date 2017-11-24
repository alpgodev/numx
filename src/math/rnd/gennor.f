c=======================================================================
c
c     GENNOR                                                 
c 
c     General Normal distribution generator N(mean,std).
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION gennor(av, sd)
c-----------------------------------------------------------------------
c
c     INPUT :
c          av  : mean of the normal distribution                  double
c          std : std of the normal distribution (>0)              double
c
c     OUTPUT :
c          gennor : generated normal deviate                      double
c
c-----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION av,sd
c
c     external functions
      DOUBLE PRECISION snorm
      EXTERNAL snorm
c     
c     executable statements
      gennor = sd*snorm() + av
      RETURN
      END
