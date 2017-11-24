c=======================================================================
c
c     subroutine MEDIAN
c
c     Median of x=[x(1),...,x(n)]
c
c     median = 
c       [x(i/2-1) + x(p/2)]/2, if p is even       
c       x(int(p/2)),           if p is odd    
c
c-----------------------------------------------------------------------
      SUBROUTINE MEDIAN ( n, x, y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of values ( n > 1)                         integer
c       x    : values (n)                                         double 
c    
c     OUTPUT 
c       y    : median                                             double
c       info : diagnostic argument                               integer
c    
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n, info
      DOUBLE PRECISION y, x(*)
c
c     
c     intrinsic functions
      INTRINSIC mod, ceiling
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      y    = 0.0
c
c     test n
      IF (n .LE. 0) THEN
        info = -2
        RETURN
      ENDIF
c
c     case n = 1
      IF (n .EQ. 1) THEN      
        y = x(1)
        RETURN
      ENDIF  
c
c     calculate median
      IF (mod(n,2) .NE. 0) THEN
        y = x(ceiling(n/2.0))
      ELSE
        y = (x(n/2)+x(n/2 + 1))/2.0
      ENDIF
c
      RETURN
      END
