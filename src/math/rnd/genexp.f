c=======================================================================
c
c     GENEXP                                                
c 
c     Exponential distribution generator
c
c     Generates a single random deviate from an exponential
c     distribution with mean AV.
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION genexp(av)
c-----------------------------------------------------------------------
c
c     INPUT :
c      AV : the mean of the exponential distribution from which
c           a random deviate is to be generated (AV >= 0)         double
c     
c     Method
c
c     Renames SEXPO from TOMS as slightly modified by BWB to use RANF
c     instead of SUNIF.
c
c     For details see:
c               Ahrens, J.H. and Dieter, U.
c               Computer Methods for Sampling From the
c               Exponential and Normal Distributions.
c               Comm. ACM, 15,10 (Oct. 1972), 873 - 882.
c
c----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION av
c     
c     external functions
      DOUBLE PRECISION sexpo
      EXTERNAL sexpo
c     
c     executable statements
 10   genexp = sexpo()*av
      RETURN
      END
