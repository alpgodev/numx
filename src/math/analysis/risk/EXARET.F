c=======================================================================
c
c     subroutine EXARET                                     
c
c     Ex-ante return := rho * w
c
c-----------------------------------------------------------------------
      SUBROUTINE EXARET ( nasset, rmean, weight, exantr, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nasset : number of asset(s)                         integer
c            rmean  : asset(s) mean returns (nasset)              double
c            weight : weights (nasset)                            double
c
c     OUTPUT 
c            exantr : ex-ante return                              double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            XV     : computing scalar product of two vectors
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nasset, info
      DOUBLE PRECISION exantr
      DOUBLE PRECISION rmean(*), weight(*)
c
c     external subroutines
      EXTERNAL XV
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     computing ex-ante return = w'*rmean
      CALL XV ( nasset, weight, rmean, exantr )
c
      RETURN
      END

