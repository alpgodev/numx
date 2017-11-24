c=======================================================================
c
c     subroutine MODVAR                                      
c
c     Cornish-Fisher (or modified) Value-at-Risk utility
c
c---------------------------------------------------------------------
      SUBROUTINE MODVAR(mean, var, skew, kurto, zc, mvar, info)
c---------------------------------------------------------------------
c
c     INPUT :
c            mean  : mean of x                                  double
c            var   : variance of x                              double
c            skew  : skewness of x                              double
c            kurto : kurtosis of x                              double
c            zc    : quantile of the normal distribution        double  
c            
c     OUTPUT :
c            mvar  : modified Value-at-Risk                     double
c            info  : diagnostic argument                       integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o 
      INTEGER info
      DOUBLE PRECISION mean, var, skew, kurto, zc, mvar
c
c     local variables
      DOUBLE PRECISION tmp, myzero
      PARAMETER (myzero = 1.E-30)
c     
c     intrinsic functions
      INTRINSIC dsqrt
c
c----------------------------------------------------------------------
c      
c     initialization
      info = 0
      mvar = 0.0    
c      
c     test if variance is too small
      IF (var .LT. myzero) THEN
         info = -104
         RETURN
      ENDIF
c          
c     executable statement  
      tmp = zc*zc 
      mvar = dsqrt(var)* (- zc + (1./6)*(tmp-1)*skew  - 
     &     (1./24)*(tmp*zc - 3*zc)*kurto  + 
     &     (1 ./36)*(2*tmp*zc - 5*zc)*skew*skew) + mean
c      
      RETURN
      END
  

      
