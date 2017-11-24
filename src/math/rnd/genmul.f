c=======================================================================
c
c     GENMUL                                                
c 
c     Generate an observation from the multinomial distribution
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      SUBROUTINE genmul(n,p,ncat,ix)
c-----------------------------------------------------------------------      
c
c     INPUT :
c      N : number of events that will be classified into one of
c           the categories 1..NCAT                               integer
c      P : Vector of probabilities.  P(i) is the probability that
c           an event will be classified into category i.  Thus, P(i)
c           must be [0,1]. Only the first NCAT-1 P(i) must be defined
c           since P(NCAT) is 1.0 minus the sum of the first
c           NCAT-1 P(i).                                 DOUBLE (NCAT-1)
c      NCAT : number of categories.  Length of P and IX           integer
c
c     OUTPUT :
c      IX : Observation from multinomial distribution.  All IX(i)
c            will be nonnegative and their sum will be N.   INTEGER (NCAT)
c
c     Method
c
c     Algorithm from page 559 of Devroye, Luc 1986
c     Non-Uniform Random Variate Generation. Springer-Verlag, New York
c
c------------------------------------------------------------------------
c
c     scalar arguments
      INTEGER n,ncat
c     
c     array arguments
      DOUBLE PRECISION p(*)
      INTEGER ix(*)
c     
c     local scalars
      DOUBLE PRECISION prob,sum
      INTEGER i,icat,ntot
c     
c     external functions
      INTEGER ignbin
      EXTERNAL ignbin
c     
c     intrinsic functions
      INTRINSIC abs
c     
c     executable statements
c     check arguments
c     see Rand.c 
c     initialize variables
      ntot = n
      sum = 1.0
      DO 20,i = 1,ncat
          ix(i) = 0
   20 CONTINUE
c
c     generate the observation
      DO 30,icat = 1,ncat - 1
          prob = p(icat)/sum
          ix(icat) = ignbin(ntot,prob)
          ntot = ntot - ix(icat)
          IF (ntot.LE.0) RETURN
          sum = sum - p(icat)
   30 CONTINUE
      ix(ncat) = ntot
c
      RETURN
      END
