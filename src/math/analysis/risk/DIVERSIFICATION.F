c=======================================================================
c     subroutine DIVERSIFICATION                                   
c
c     % of Diversifiable risk eliminated 
c     
c     % risk eliminated = 0.5*(1 - ipc)
c
c     where ipc is the intra-portfolio correlation
c 
c----------------------------------------------------------------------
      SUBROUTINE DIVERSIFICATION ( n, corr, w, coef, info)
c----------------------------------------------------------------------
c
c     INPUT 
c       n     : portfolio size (n>=1)                          integer
c       corr  : correlation matrix (n*n)                        double
c       w     : weights (n)                                     double
c
c     OUTPUT 
c       coef  : intra-portfolio correlation coefficient         double
c       info  : diagnostic argument                            integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION coef, w(*), corr(n,*)
c
c     local variables
      DOUBLE PRECISION ipc
c
c-----------------------------------------------------------------------
c
c     initialization
      coef = 0.0
      info = 0
c
c----------------------------------------------------------------------
c
c     intra-portfolio correlation
      CALL EXAIPC ( n, corr, w, ipc, info)
      IF (info .LT. 0) RETURN
c
c     % of risk eliminated      
      coef = 0.5*(1.0 - ipc)     
c
      RETURN
      END
