c=======================================================================
c
c     subroutine APTVOL                                      
c
c     Robusts alpha and beta (APT method) vectors of an universe for 
c     a given set of factors and volatility of regression error
c
c-----------------------------------------------------------------------
      SUBROUTINE aptvol ( n, p, q, Y, X, iwork, dwork,
     &                    alpha, beta, vol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of value(s) (n >= 1)                integer
c            p      : number of asset(s)                         integer
c            q      : number of factors                          integer
c            Y      : asset(s) values (n*p)                       double
c            X      : factor(s) values (n*q)                      double
c
c     WORKSPACE 
c            iwork  : q + 1                                      integer
c            dwork  : (q + 1)*( n + 2*p + q + 2 ) + n*p           double
c
c     OUTPUT 
c            alpha  : alpha coefficient(s) (p)                    double
c            beta   : beta coefficent(s) (q*p)                    double
c            vol    : std of errors (p)                           double
c            info   : diagnostic argument                        integer
c
c     CALL   
c            CALAPT : robusts APT method
c            RHAPT  : historical assets values as a function
c                     of values of factors with the APT regression
c            DM     :
c            HVOLAT :
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, q, info
      DOUBLE PRECISION X(*), Y(*), alpha(*), beta(q,*), vol(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piw, pdw, pdrs
c
c     external subroutines
      EXTERNAL calapt, RHAPT, DM, HVOLAT
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for CALAPT workspace who needs (q+1)
c
c     Total size of iwork array = q + 1
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdrs = 1
c     pdrs : pointer for matrix of simulated returns, so ( n*p ) more
      pdw = pdrs + ( n*p )
c     pdw  : pointer for CALAPT workspace who needs
c            (q + 1)*(n + 2*p + q + 2) 
c
c     Total size of dwork array = (q + 1)*( n + 2*p + q + 2 ) + n*p
c
c-----------------------------------------------------------------------
c
c     construction of matrix F
      CALL calapt ( n, p, q, Y, X, iwork(piw), dwork(pdw),
     &              alpha, beta, info )
      IF (info .lt. 0) RETURN
c     
c     simulated history of returns and alpha/beta
      CALL RHAPT ( n, p, q, X, alpha, beta, dwork(pdrs) )
c    
c     regression error
      CALL DM ( n, p, Y, dwork(pdrs), dwork(pdrs) )
c     
c     regression std error (for each asset)
      CALL HVOLAT ( n, p, dwork(pdrs), vol, info )
c
      RETURN
      END

