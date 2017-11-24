c=======================================================================
c
c     subroutine EXARVA                                     
c
c     Ex-ante Variance
c
c-----------------------------------------------------------------------
      SUBROUTINE EXARVA ( n, cov, w, var, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c        n      : number of asset(s)                             integer
c        cov    : covariance matrix  (n*n)                        double
c        w      : positions/weights (n)                           double
c
c     OUTPUT 
c        var    : variance                                        double
c        info   : diagnostic argument                            integer
c
c     CALL   
c        OVTMCV : M'*D*M = scalar
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION var, w(*), cov(*) 
c
c-----------------------------------------------------------------------
c
c     initialization
      var  = 0.0 
      info = 0
c
c     variance = w' * Cov * w
      CALL OVTMCV (n, cov, w, var )      
c
      RETURN
      END
