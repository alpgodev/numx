c=======================================================================
c     
c     subroutine SIMSTRACE
c
c     Subroutine used by BFGS to compute the value of the dual function
c     and its gradient specific for Trace constraint
c
c-----------------------------------------------------------------------
      SUBROUTINE simstrace ( indic, simext, mct, ydual, funct, grad,
     &                       iwork, dwork)
c-----------------------------------------------------------------------
c
c     INPUT 
c       indic  : = 4, to compute function and gradients          integer
c       simext : entry point of an external subroutine
c                provided by the user (not used here)
c       mct    : number of variables                             integer
c       ydual  : dual solution (mct)                              double
c
c     OUTPUT 
c       funct  : dual function value                              double
c       grad   : gradient value mct)                              double
c
c     WORKSPACE 
c       iwork  : 12*n + 2                                        integer
c                contents in output :
c                n (size of the problem)    in iwork(1)
c                info (errors of SOLSTRACE) in iwork(2)
c       dwork  : n*(6*n + 27)                                     double
c                contents in input :
c                matrix C(n*n)
c     CALL   
c       SOLSTRACE : computing the solution of SDLS optimization from the
c                   dual solution specific to Trace constraint
c       NMS       : computing the Frobenius norm of a vectorized symmetric matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, mct
      DOUBLE PRECISION funct, ydual, grad
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER n, pn, pc, px, pdw, piw, pinfo
      DOUBLE PRECISION normX, traceX, traceC, yb
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pn     = 1
c     pn     : pointer for n the size of the problem, so 1 more
      pinfo  = pn + 1
c     pinfo  : pointer for info (errors of SOLSTRACE), so 1 more
      piw    = pinfo + 1
c     piw    : pointer for workspace used by SOLSTRACE, so 12*n more
c
c     Total size of iwork array (12*n + 2)
c
c     initializations
      n   = iwork(pn)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pc = 1
c     pc    : pointer for the matrix C (n*n)
      px  = pc + ( n*n )
c     px    : pointer for the solution matrix X ( n*n ) 
      pdw = px + ( n*n )
c     pdw   : pointer for workspace used by SOLSTRACE (n*(4*n + 27)
c
c     Total size of dwork array = n*(6*n + 27)
c
c-----------------------------------------------------------------------
c
c     projection X = pK( C + A*y )
      CALL solstrace ( n, dwork(pc), ydual, iwork(piw), dwork(pdw),
     &                 dwork(px), iwork(pinfo) )
c
c     computing ||x||**2
      CALL NM ( n, n, dwork(px), normX )
c
c     computing scalar product ydual by b (Trace(C))
      CALL TM ( n, dwork(pc), traceC )
      yb = traceC*ydual
c
c     computing function value
c     f = -y'b + 0.5*||pK(C+A*y)||**2
c       = -y'b + 0.5*||X**2||
      funct = -yb + 0.5*(normX*normX)
c
c     computing gradient, grad = A*pK( C+A*y ) - b
c                              = Trace(X) - Trace(C)
c     computing trace de X
      CALL TM ( n, dwork(px), traceX )
      grad = traceX - traceC
c
      n=indic+mct
      RETURN
      END
