c=======================================================================
c
c     subroutine EXPVRI                                      
c
c     Upper bound of the Ex-post Value-at-Risk
c
c------------------------------------------------------------------------
      SUBROUTINE EXPVRI ( n, x, alpha, valris, info)
c------------------------------------------------------------------------
c
c     INPUT 
c        n     : number of point(s)                               integer
c        x     : value(s) (n)                                      double
c        alpha : confidence probability level (0<alpha<1)          double
c
c     OUTPUT 
c        valris : upper bound Value-at-Risk                        double
c        info   : = 0 successful exit                             integer
c
c     CALL   
c        EXPRET   : ex-post return
c        EXPVOL   : ex-post volatility
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n, info
      DOUBLE PRECISION alpha, valris
      DOUBLE PRECISION x(*)
c
c     local variables
      DOUBLE PRECISION vol, ret, eps
      PARAMETER ( eps = 1.E-8 )
c
c     external subroutines
      EXTERNAL EXPRET, EXPVOL      
c     
c     intrinsic functions
      INTRINSIC sqrt, min, max  
c
c-----------------------------------------------------------------------
c
c     initialization
      valris = 0.0
      info   = 0
c
c     test confidence level:  eps <= alpha <= 1-eps
      IF ( (alpha .LT. eps) .or. (alpha .GT. (1. - eps)) ) THEN
         info = -102
         RETURN
      ENDIF
c
c     ex-post return
      CALL EXPRET ( n, x, ret, info )
      IF (info .LT. 0) RETURN
c
c     ex-post volatility
      CALL VOLAT ( n, x, vol, info )
      IF (info .LT. 0) RETURN
c
c     computing Value-at-Risk upper bound    
      valris = - vol * sqrt( alpha/(1.-alpha) ) + ret
c      valris = max( min(valris, 1.), -1. )
c
      RETURN
      END
c
c=======================================================================
