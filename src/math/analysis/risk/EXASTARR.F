c=======================================================================
c
c     subroutine EXASTARR                                      
c
c     Ex-ante Conditionnal Sharpe ratio (STARR)
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
c     Also known as STARR ratio, CVaR Ratio.
c
c-----------------------------------------------------------------------
      SUBROUTINE EXASTARR(n, cov, rho, w, rfr, alpha, starr, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n        : number of asset(s)                            integer
c       cov      : covariance matrix (n*n)                        double
c       rho      : expected returns (n)                           double
c       w        : weights (n)                                    double 
c       rfr      : risk-free rate                                 double
c       alpha    : confidence level                               double
c
c     OUTPUT 
c       starr    : ex-ante STARR ratio                            double
c       info     : diagnostic argument                           integer
c
c     CALL   
c        EXARET  : ex-ante return
c        EXACVAR : ex-ante Conditional Value-at-Risk
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION alpha, rfr, starr, w(*), rho(*), cov(*)
c
c     local variables
      DOUBLE PRECISION mean, cvar, EPS
      PARAMETER (EPS = 1.E-8)
c
c     external subroutines
      EXTERNAL EXARET, EXACVAR  
c
c-----------------------------------------------------------------------
c
c     initialisation 
      info  = 0
      starr = 0.0
c
c     test confidence level: alpha
      IF ((alpha .LT. EPS) .OR. (alpha .GT. (1. - EPS))) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-ante return
      CALL EXARET ( n, rho, w, mean, info )
      IF (info .LT. 0) RETURN
c
c     ex-ante Conditional Value-at-Risk
      CALL EXACVAR ( n, cov, rho, w, alpha, cvar, info)
      IF (info .LT. 0) RETURN
c
c     STARR ratio
      starr = (mean - rfr)/ABS(cvar)
      RETURN
      END
