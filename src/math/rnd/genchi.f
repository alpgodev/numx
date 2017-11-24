c=======================================================================
c
c     GENCHI                                                 
c 
c     Chi-square distribution generator
c
c     Generates random deviate from the distribution of a chisquare
C     with DF degrees of freedom random variable.
c
c-----------------------------------------------------------------------
c
c     Copyright (c) 2014 NumX
c     All rights reserved.
c 
c     This software is the confidential and proprietary information
c     of NumX. You shall not disclose such Confidential
C     Information and shall use it only in accordance with the terms
c     of the licence agreement you entered into with NumX.
c
c     author: Yann Vernaz
c
c----------------------------------------------------------------------- 
      DOUBLE PRECISION FUNCTION genchi(df)
c-----------------------------------------------------------------------
c
c     INPUT :
c      DF : degrees of freedom of the chisquare (DF>0)           double
c
c     Method - Uses relation between chisquare and gamma.
c
c----------------------------------------------------------------------
c
c     scalar arguments 
      DOUBLE PRECISION df
c     
c     external functions
c      DOUBLE PRECISION gengam
c      EXTERNAL gengam
      DOUBLE PRECISION sgamma
      EXTERNAL sgamma
c     
c     executable statements
c     changed this to call sgamma directly
c   10 genchi = 2.0*gengam(1.0,df/2.0)
 10   genchi = 2.0*sgamma(df/2.0)
      RETURN
      END

