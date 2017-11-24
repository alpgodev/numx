c=======================================================================
c
c     function ranf
c
c     This function Generates Uniform[0,1] random numbers. 
c     The function has no arguments. 
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
c
      FUNCTION ranf()
c
      IMPLICIT NONE
      REAL ranf
c
c-----------------------------------------------------------------------
c      
c     Pseudorandom numbers from the uniform distribution 
c     within the range 0 <= x <=1
#IFDEF INTELFOR 
      CALL init_random_seed()
      CALL RANDOM_NUMBER(ranf)
#ELSE 
      ranf = rand()      
#ENDIF        
      RETURN
      END
