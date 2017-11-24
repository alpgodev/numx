c=======================================================================
c
c     GENNCH                                                
c 
c     Generate random value of Noncentral chi-square variable
c
c     Generates random deviate from the distribution of a noncentral
c     chisquare with DF degrees of freedom and noncentrality parameter
c     XNONC.
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
      DOUBLE PRECISION FUNCTION gennch(df, xnonc)
c-----------------------------------------------------------------------
c
c     INPUT :
c      df    : degrees of freedom of the chisquare (>=1.0)       double
c      xnonc : noncentrality parameter of the chisquare (>=0)    double
c
c      Method
c
c     Uses fact that  noncentral chisquare  is  the  sum of a  chisquare
c     deviate with DF-1  degrees of freedom plus the  square of a normal
c     deviate with mean sqrt(XNONC) and standard deviation 1.
c
c-----------------------------------------------------------------------
c
c     scalar arguments
      DOUBLE PRECISION df, xnonc
c     
c     external functions 
c     changed these to call GENCHI and GENNOR directly
c      DOUBLE PRECISION genchi, gennor
c      EXTERNAL genchi, gennor
      DOUBLE PRECISION sgamma, snorm
      EXTERNAL sgamma, snorm
c     
c     intrinsic functions
      INTRINSIC sqrt
c     
c     changed abort to df < 1, and added case: df = 1 
c     executable statements
c     changed this to call GENCHI and GENNOR directly
c      gennch = genchi(df-1.0) + gennor(sqrt(xnonc),1.0)**2
 10   IF (df.GE.1.000001) GO TO 20
c
c     case DF = 1.0
      gennch = (snorm() + sqrt(xnonc))**2
      GO TO 30
c      
c     case DF > 1.0
 20   gennch = 2.0*sgamma((df-1.0)/2.0) + (snorm() + sqrt(xnonc))**2
 30   RETURN     
      END
c
c=======================================================================