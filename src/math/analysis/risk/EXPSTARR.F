c=======================================================================
c
c     subroutine EXPSTARR                                     
c
c     Ex-ante Conditionnal Sharpe ratio, 
c     also known as STARR ratio, CVaR Ratio.
c     
c     CSR = {r - rfr} / CVaR
c
c     r    : portfolio return
c     rfr  : riskfree rate
c     CVaR : Conditional VaR, CVaR = E[-r | r <= -VaR]
c
c     Assuming a normal distribution, CVaR can be calculated as follows:
c     CVaR = avr + vol * d(N(1-alpha))/ alpha
c
c     avr        : average return
c     vol        : volatility
c     d(.)       : density function of the standard normal distribution
c     alpha      : confidence level
c     N(1-alpha) : 1-alpha quantil of the standard normal distribution
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPSTARR ( n, x, rfr, alpha, starr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n     : number of return(s) (n > 1)                     integer
c       x     : return(s) (n)                                    double
c       rfr   : risk-free rate                                   double
c       alpha : probability confidence level (0 < alpha < 1)     double
c
c     OUTPUT 
c       starr : ex-post STARR ratio                              double
c       info  : diagnostic argument                             integer
c
c     CALL   
c       EXPRET  : ex-post mean return
c       VOLAT   : ex-post volatility
c       dinvnr  : cumulative normal distribution inverse
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION rfr, alpha, starr
      DOUBLE PRECISION x(*)
c
c     local variables      
      DOUBLE PRECISION meanx, cvar, EPS
      PARAMETER ( EPS=1.E-8 )
c
c     external subroutines
      EXTERNAL EXPRET, EXPCVAR
c
c-----------------------------------------------------------------------
c
c     initialization
      starr = 0.0
      info = 0
c
c     test confidence level
      IF ((alpha .LT. EPS) .OR. (alpha .GT. (1.-EPS))) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-post mean return
      CALL EXPRET ( n, x, meanx, info )
      IF (info .LT. 0) RETURN
c
c     ex-post Conditional Value-at-Risk
      CALL EXPCVAR ( n, x, alpha, cvar, info)
c
c     STARR ratio
      starr = (meanx-rfr)/ABS(cvar)
      RETURN
      END
