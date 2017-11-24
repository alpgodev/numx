c=======================================================================
c
c     subroutine  SDLSTRACE                                    
c
c     Semi-Definite Least Square optimization
c     
c     min || X - C ||**2
c      x 
c     s.t. 
c        Trace(X) = Trace(C)
c        X >= alpha*Id         with 0.0 =< alpha < 1.0  
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlstrace ( n, C, epsbfg, alpha, iwork, dwork, X, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : matrix size                                    integer
c       C      : matrix (n*n) to optimize                        double
c       epsbfg : BFGS precision stop test                        double
c       alpha  : 0.0 =<alpha < 1.0                               double
c
c     WORKSPACE 
c       iwork  : 12*n + 5                                       integer
c       dwork  : n*(7*n + 27) + 6                                double       
c
c     OUTPUT 
c       X      : optimal matrix (n*n)                            double
c       info   : diagnostic argument                            integer
c
c     CALL   
c        SEMDX    : computing the sum of the diagonal elements of a
c                   matrix with a scalar : B = A + x*I
c                  ( A square matrix(n*n), x scalar,
c                    gives B square matrix(n*n) )
c        PMX      : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YM       : copy a vectorized matrix in a vectorized matrix
c        BFGSSDLS : quasi-Newton optimizer with BFGS method
c                   ( short call version for SDLS )
c        SIMSTRACE: subroutine used by BFGS to compute the value of the
c                   dual function and its gradient specific to trace pb
c        SOLSTRACE: computing the solution of SDLS optimization from the
c                   dual solution   x = pK( C + Ai.y(i) ) 
c                   specific to obtain a correlation matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simstrace, gesterr
c
c     i/o arguments
      INTEGER n, info
      DOUBLE PRECISION epsbfg, alpha, C(*), X(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mct, totitr, totsim, piw, pinfo, pdw, pdcv, pdx, pn
      DOUBLE PRECISION ydual, binf, bsup, funct, grad, eps, infini, 
     &                 traceC
      PARAMETER ( eps = 1.e-15, infini = 1.d20 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0   
      mct  = 1
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw   = 1
c     piw   : pointer for local workspaces,
c             BFGSSDLS needs (2*mct+1)
c                      + communication with SIMSTRACE (12*n+2),
c                      so ((2*mct+1)+(12*n+2))
c             SOLSTRACE needs (12*n)
c             so, union of spaces : ((2*mct+1)+(12*n+2))
c             for trace constraint mct = 1,
c             so : (12*n+5)
c
c     Total size of iwork array = 12*n + 5
c
      pn     = piw + (2*mct+1)
c     pn     : communication pointer with SIMSTRACE for n
      pinfo  = pn + 1
c     pinfo  : communication pointer with SIMSTRACE for info
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdx    = 1
c     pdxs   : pointer for matrix X (n*n)
      pdw    = pdx + ( n*n )
c     pdw    : pointer for local workspaces,
c             BFGSSDLS needs ( mct*(mct+11)/2 ) -> 6
c                      + communication with SIMSTRACE n*(6*n + 27)
c             SOLSTRACE needs n*(4*n + 27)
c             so, union of spaces: n*(6*n + 27) + 6
c
c     Total size of dwork array = n*(7*n + 27) + 6
c 
      pdcv  = pdw +  6
c     pdcv  :  communication pointer with SIMSTRACE
c              for the vectorized symmetric matrix C(n*n),
c
c-----------------------------------------------------------------------
c
c     problem feasibility: n*alpha <= Trace(C) 
      CALL TM ( n, C, traceC )
      IF (traceC .LT. (n*alpha)) THEN
        info = -1300
        RETURN
      ENDIF
c      
c     construction of dwork vector (communication with simstrace)
c     chgt. variable: C = C - alpha*I (alpha > 0)
      IF (alpha .GT. eps) THEN
         CALL SEMDX ( n, C, -alpha, dwork(pdcv) )
      ELSE
         CALL YM ( n, n, C, dwork(pdcv) )
      ENDIF
c
c     initialization of ydual
      ydual = 0.
c
c     initializations of inf/sup bounds
      binf = -infini
      bsup =  infini
c
c     construction of iwork vector (communication with simstrace)
      iwork(pn)    = n
      iwork(pinfo) = 0
c
c     non-linear optimization (BFGS) -> dual solution
c
c     (the second parameter is not used here, we can replace GESTERR
c      by an other name of subroutine)
      CALL bfgssdls ( simstrace, gesterr, mct,ydual, epsbfg, binf, bsup,
     &                iwork(piw), dwork(pdw),
     &                funct, grad, totitr, totsim, info)
c
c     gestion of SIMSTRACE errors
      info = iwork(pinfo)
c
c     construction of the matrix solution from dual solution
c     projection on the SDP cone
      CALL solstrace ( n, dwork(pdcv), ydual,
     &                 iwork(piw), dwork(pdw), dwork(pdx), info )
c
c     chgt. variable: X = Z + alpha*I (alpha > 0)
      IF (alpha .GT. eps) THEN
         CALL SEMDX ( n, dwork(pdx), alpha, X )
      ELSE
         CALL YM ( n, n, dwork(pdx), X )
      ENDIF
c
c     test if Trace(C) < 0 (then priority to constraint X >= alpha*Id)
      IF (traceC .LT. eps) info = 1301
c
      RETURN
      END
