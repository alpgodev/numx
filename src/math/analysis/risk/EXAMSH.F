c=======================================================================
c
c     subroutine EXAMSH                                      
c
c     Ex-ante Cornish-Fisher Sharpe ratio (modified Sharpe ratio)
c
c-----------------------------------------------------------------------
      SUBROUTINE EXAMSH(n, p, w, x, cov, riskfr, alpha, dwork, 
     &                  examsha, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values                           integer
c            p      : number of asset(s)                         integer
c            w      : weight(s) (p)                               double 
c            x      : asset(s) values (n*p)                       double
c            cov   : covariance matrix (p*p)                      double
c            riskfr : risk-free rate                              double
c            alpha  : confidence level                            double
c
c     WORKSPACE 
c            dwork  : 3*p*( p + n + 1 )                           double
c
c     OUTPUT 
c            examsha : ex-ante Cornish-Fisher Sharpe ratio        double
c            info    : diagnostic argument                       integer
c
c     CALL   
c        EXARET  : ex-ante return
c        EXAMVA  : ex-ante Cornish-Fisher Value-at-Risk
c        MODSHA  : modified Sharpe ratio
c
c-----------------------------------------------------------------------
c      
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION alpha, riskfr, examsha
      DOUBLE PRECISION w(*), x(*), cov(*)
c
c     workspaces      
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdwork, pdrho, pdcov, pdwcov, pdrvar
      DOUBLE PRECISION examvar, exarho, EPS
      PARAMETER (EPS = 1.E-8)
c
c     external subroutines
      EXTERNAL COVMVM, EXARET, EXAMVA, MODSHA  
c
c-----------------------------------------------------------------------
c
c     initialisation 
      examsha = 0.0
c      
c     pointers for double work space : dwork
c     -------------------------------------- 
      pdwork = 1
      pdrho  = pdwork 
c     pdrho : pointer for average return, so ( p ) more
      pdcov  = pdrho + ( p )
c     pdcov : pointer for cov matrix, so ( p*p ) more     
      pdwcov = pdcov + ( p*p)
c     pdwcov : pointer dwork COVMVM, so ( p*n ) more       
      pdrvar = pdwcov + ( p*n )
c     pdrvar : pointer dwork EXAMVA p * ( 6*p + n + 2 ) 
c
c     Total dwork = p*( 7*p + 2*n + 3 )  
c
c-----------------------------------------------------------------------
c
c     test confidence level: alpha
      IF ((alpha .LT. EPS) .or. (alpha .GT. (1. - EPS))) THEN
         info = -102
         RETURN
      ENDIF
c
c     covariance matrix and average return 
      CALL COVMVM ( n, p, x, dwork(pdwcov), dwork(pdrho), 
     &            dwork(pdcov), info)
      IF (info .LT. 0) RETURN
c
c     ex-ante return
      CALL EXARET ( p, dwork(pdrho), w, exarho, info )
      IF (info .LT. 0) RETURN
c
c     Cornish-Fisher (modified) ex-ante Value-at-Risk
      CALL EXAMVA(n, p, w, x, cov, alpha, dwork(pdrvar), examvar, info)
      IF (info .LT. 0) RETURN
c
c     Cornish-Fisher (modified) Sharpe ratio
      CALL MODSHA(exarho, riskfr, examvar, examsha, info)
c
      RETURN
      END
