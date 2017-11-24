c=======================================================================
c
c     subroutine GREPS
c
c     Computing groups in a sorted vector for a sensibility
c
c-----------------------------------------------------------------------
      SUBROUTINE greps ( N, X, eps, NG, XG, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        N    : size of the vector                               integer
c        X(N) : input vector (N), in decrease order               double
c        eps  : groups sensibility                                double
c
c     OUTPUT 
c        NG   : number of groups                                 integer
c        XG   : groups vector (NG), max = N                      integer
c        info : diagnostic argument                              integer
c
c     METHOD 
c            for an example of vector ( N = 7 ) :   
c            X = [ 5.2, 5.1, 3.5, 3.4, 3.3, 1.7 1.3 ]    
c            and eps = 0.15, the routine computes    
c            NG = 4 and XG = [ 2, 3, 1, 1 ]    
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o parameters
      INTEGER N, NG, info, XG(*)
      DOUBLE PRECISION eps, X(*)
c
c     local variables
      INTEGER i
c
c     intrinsic functions
      INTRINSIC ABS
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     one group by default 
      NG = 1
      XG(NG) = 1
c
c     loop      
      DO i = 2,N
         IF ( ABS( X(i) - X(i - 1) ) .LT. eps ) THEN
            XG(NG) = XG(NG) + 1
         ELSE
            NG = NG + 1
            XG(NG) = 1
         ENDIF
      ENDDO
      RETURN
      END
