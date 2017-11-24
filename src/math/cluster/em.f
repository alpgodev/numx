c=======================================================================
c
c     subroutine EM                                    
c
c     EM Algorithm (Cluster Analysis) 
c
c-----------------------------------------------------------------------
      SUBROUTINE em ( itermax, n, d, k, x, iwork, dwork,
     &                pk, nk, muk, covk, label, proba, f, like, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       itermax  : number of iterations max.                     integer
c       n        : number of points (n>1)                        integer
c       d        : dimension (d>0)                               integer
c       k        : number of class (k>0)                         integer
c       x        : data (n*d)                                     double
c
c     WORKSPACE
c       iwork    : d                                             integer
c       dwork    : d*(3*d + 9)                                    double
c
c     INPUT/OUTPUT 
c       pk       : proportion of each group (k)                   double
c       nk       : number of each group (k)                      integer
c       muk      : mean of each group (d*k)                       double
c       covk     : cov. matrix of each group (d)*(d*k)            double
c       label    : partition labels (n*k)                        integer
c       proba    : conditionale probability (n*k)                 double
c       f        : probability (n*k)                              double
c       like     : likelihood value                               double                                      
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER itermax, n, d, k, nk(*), label(n,*), info
      DOUBLE PRECISION x(n,*), pk(*), proba(n,*), muk(d,*), covk(d,*), 
     &                 f(n,*), like
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER l, iter, piw, pdw
c
c     external subroutines
      EXTERNAL estep, mstep
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw = 1
c     piw : pointer for ESTEP, so ( d ) more
c
c     Total size of dwork array = ( d ) 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdw = 1
c     pdw : pointer for ESTEP, so ( d*(3*d + 9) ) more
c
c     Total size of dwork array = d*(3*d + 9) 
c
c------------------------------------------------------------------------
c
c     loop EM
      DO iter = 1, itermax
c
c       E step (expectation)
        CALL estep ( n, d, k, x, pk, muk, covk,
     &               iwork(piw), dwork(pdw), label, proba, f, like,info)
        IF (info .LT. 0) RETURN
c
c       M step (maximization)
        CALL mstep ( n, d, k, x, label, proba,
     &               pk, nk, muk, covk, info )
        IF (info .LT. 0) RETURN
        DO l = 1,k
            IF (nk(l) .EQ. 0) THEN
                like = -1.E+15
                RETURN
            ENDIF
        ENDDO
      ENDDO
c
      RETURN
      END
     
