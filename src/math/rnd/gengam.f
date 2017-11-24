c=======================================================================
c
c     GENGAM                                                 
c 
c     Gamma distribution generator
c
c     Generates random deviates from the Gamma distribution whose
c     density is
c
c          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION gengam(a, r)
c-----------------------------------------------------------------------      
c
c     INPUT :
c      a : location parameter of Gamma distribution (a>0)         double
c      r : shape parameter of Gamma distribution (r>0)            double
c
c     Method - for details see:
c               (Case R >= 1.0)
c               Ahrens, J.H. and Dieter, U.
c               Generating Gamma Variates by a
c               Modified Rejection Technique.
c               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
c     Algorithm GD
c
c               (Case 0.0 < R < 1.0)
c               Ahrens, J.H. and Dieter, U.
c               Computer Methods for Sampling from Gamma,
c               Beta, Poisson and Binomial Distributions.
c               Computing, 12 (1974), 223-246/
c     Adapted algorithm GS.
c
c----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION a, r
c     
c     external functions
      DOUBLE PRECISION sgamma
      EXTERNAL sgamma
c     
c     executable statements
 10   gengam = sgamma(r)/a
c      gengam = gengam/a
      RETURN
      END
