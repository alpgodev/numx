c=======================================================================
c
c     subroutine MODSHA                                      
c
c     Cornish-Fisher (or modified) Sharpe ratio utility
c
c----------------------------------------------------------------------
      SUBROUTINE MODSHA(rho, riskfr, modvar, modsh, info)
c----------------------------------------------------------------------
c
c     INPUT :
c            rho    : mean return                                double
c            riskfr : risk-free rate                             double
c            modvar : modified Value-at-Risk                     double  
c            
c
c     OUTPUT :
c            modsh : modified Sharpe ratio                       double
c            info  : diagnostic argument                        integer
c     
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER info
      DOUBLE PRECISION rho, riskfr, modvar, modsh
c
c     local variables      
      DOUBLE PRECISION myzero
      PARAMETER (myzero = 1.E-30)
c     
c     intrinsic functions
      INTRINSIC abs       
c
c----------------------------------------------------------------------
c      
c     initialization
      info  = 0
      modsh = 0.0
c
c     test if modified Value-at-Risk is too small 
      IF (ABS(modvar) .lt. myzero) THEN
         info = -1
         RETURN
      ENDIF
c     
c     executable statement
      modsh = (rho-riskfr)/ABS(modvar)
c
      RETURN
      END
