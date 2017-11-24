c=======================================================================
c
c     subroutine EXPCVAR                                    
c
c     Ex-post normal (i.e. Gaussian) conditional Value-at-Risk
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPCVAR ( n, x, alpha, cvar, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n     : number of return(s) (n > 1)                 integer
c            x     : return(s) (n)                                double
c            alpha : probability confidence level (0 < alpha < 1) double
c
c     OUTPUT 
c            cvar   : Gaussian conditional Value-at-Risk          double
c            info   : diagnostic argument                        integer
c
c     CALL   
c          EXPRET  : ex-post mean return
c          VOLAT   : ex-post volatility
c          dinvnr  : cumulative normal distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha, cvar
      DOUBLE PRECISION x(*)
c
c     local variables      
      DOUBLE PRECISION zc, p, q, volx, meanx, coeff, eps, pi
      PARAMETER ( eps = 1.E-8, pi = 3.1415926535898 )
c
c     external functions
      DOUBLE PRECISION dinvnr
      EXTERNAL dinvnr
c
c     external subroutines
      EXTERNAL EXPRET, VOLAT
c     
c     intrinsic functions
      INTRINSIC sqrt, exp, min, max        
c
c-----------------------------------------------------------------------
c
c     initialization
      cvar = 0.0
      info = 0
c
c     test confidence level: alpha
      IF ((alpha .LT. eps) .OR. (alpha .GT. (1.-eps))) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-post mean return
      CALL EXPRET ( n, x, meanx, info )
      IF (info .LT. 0) RETURN
c
c     ex-post volatility
      CALL VOLAT ( n, x, volx, info )
      IF (info .LT. 0) RETURN
c
c     computing inverse of the N(0,1) cumulative distribution
c     call the function dinvnr(p,q)
      p = alpha
      q = 1. - p
      zc = dinvnr(p,q)
c
c     coeff := [sqrt(2*pi)*exp[(zc/2)^2]*(1-alpha)]^(-1)
      coeff = sqrt(2.0*pi)*exp(((zc*zc)/4.0))*(1.0 - alpha)
      coeff = 1.0/coeff
c     
c     computing Normal Conditional Value-at-Risk 
c     CVaR := max[0, min(1, -coeff*volatility - mean)] 
      cvar = - volx * coeff + meanx 
c      cvar = min(cvar, 1.)
c      cvar = max(cvar, -1.)
c
      RETURN
      END
