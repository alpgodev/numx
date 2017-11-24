c=======================================================================
c     subroutine EXACC                                      
c
c     Concentration coefficient 
c     
c                               1    
c     concentration coef. = ---------
c                           || w ||^2           
c 
c     where ||.|| is the L2-norm
c
c----------------------------------------------------------------------
      SUBROUTINE EXACC ( n, w, coef, info)
c----------------------------------------------------------------------
c
c     INPUT 
c       n     : portfolio size (n>=1)                          integer
c       w     : weights (n)                                     double
c
c     OUTPUT 
c       coef  : concentration coefficient                       double
c       info  : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION coef, w(*)
c
c     local variables
      DOUBLE PRECISION norm, EPS
      PARAMETER ( EPS = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialization
      coef = 0.0
      info = 0
c
c----------------------------------------------------------------------
c                                    _________________
c     L2-norm of a vector: ||w|| = \/ sum of w(i)*w(i)
      CALL NV ( n, w, norm )
      IF (norm .GT. EPS) THEN
        coef = 1.0/(norm*norm)
      ELSE
        info = -1
      ENDIF      
c
      RETURN
      END

