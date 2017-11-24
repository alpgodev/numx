c=======================================================================
c
c     SEXPO
c 
c     Generate Standard Exponential distribution
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
      DOUBLE PRECISION FUNCTION sexpo()
c-----------------------------------------------------------------------      
c
c     Method - for details see:                                                 
c                                                                      
c     AHRENS, J.H. AND DIETER, U.                            
c     COMPUTER METHODS FOR SAMPLING FROM THE                 
c     EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  
c     COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               
c                                                                      
c     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
c     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
c
c     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
c     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
c
c-----------------------------------------------------------------------
c
c     local scalars
      DOUBLE PRECISION a, q1, u, umin, ustar
      INTEGER i
c     
c     local arrays
      DOUBLE PRECISION q(8)
c     
c     external functions 
      DOUBLE PRECISION ranf
      EXTERNAL ranf
c     
c     equivalences
      EQUIVALENCE (q(1), q1)
c
c     save statement
      SAVE q
c    
c     data statements
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
c
   10 a = 0.0
      u = ranf()
      GO TO 30

   20 a = a + q1
   30 u = u + u
c     changed the following to reflect the true algorithm and
c     prevent unpredictable behavior if U is initially 0.5.
c      IF (u.LE.1.0) GO TO 20
      IF (u.LT.1.0) GO TO 20
   40 u = u - 1.0
      IF (u.GT.q1) GO TO 60
   50 sexpo = a + u
      RETURN
c
   60 i = 1
      ustar = ranf()
      umin = ustar
   70 ustar = ranf()
      IF (ustar.LT.umin) umin = ustar
   80 i = i + 1
      IF (u.GT.q(i)) GO TO 70
   90 sexpo = a + umin*q1
      RETURN
      END
