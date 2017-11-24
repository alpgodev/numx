c=======================================================================
c
c     subroutine EXPTER                                      
c
c     Ex-post Tracking Error
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPTER ( n, x, retben, dwork, tracke, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c          n      : size of x (n > 1)                            integer
c          x      : returns (n)                                   double
c          retben : returns of benchmark (n)                      double
c
c     WORKSPACE 
c          dwork  : vector ( n )                                  double
c
c     OUTPUT 
c          tracke : tracking error                                double
c          info   : diagnostic argument                          integer
c
c     CALL   
c          VOLAT  : volatility (standard deviation)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION tracke
      DOUBLE PRECISION x(*), retben(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i
c
c     external subroutines
      EXTERNAL VOLAT        
c
c-----------------------------------------------------------------------
c
c     initialization
      tracke = 0.0
      info   = 0
c
c     computing ( return portfolio - return benchmark )
      DO i = 1,n
         dwork(i) = x(i) - retben(i)
      ENDDO
c
c     computing volatility
      CALL VOLAT ( n, dwork, tracke, info )
c 
      RETURN
      END
c
c=======================================================================
