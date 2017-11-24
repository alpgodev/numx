c=======================================================================
c
c     function init_random_seed                                        
c
c     Initialize a pseudo-random number sequence
c
c-----------------------------------------------------------------------
c
c     Copyright (c) 2012 NumX
c     515, route de Savoie, 38114 Allemont, France
c     All rights reserved.
c 
c     This software is the confidential and proprietary information
c     of NumX. You shall not disclose such Confidential
c     Information and shall use it only in accordance with the terms
c     of the licence agreement you entered into with NumX.
c
c     Author : Yann Vernaz
c
c-----------------------------------------------------------------------
      SUBROUTINE init_random_seed()
c-----------------------------------------------------------------------
c            
      REAL, DIMENSION(100) :: rnd
      INTEGER              :: isize,idate(8)
      INTEGER,ALLOCATABLE  :: iseed(:)

      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize)
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed)
      iseed = iseed * (idate(8)-500)  ! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)

c      CALL RANDOM_NUMBER(RND)

      DEALLOCATE( iseed )
c      
      END SUBROUTINE
