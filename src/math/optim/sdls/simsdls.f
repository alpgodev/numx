c=======================================================================
c
c     SIMSDLS
c
c     Subroutine used by BFGS to compute the value of the dual function
c     and its gradient
c
c-----------------------------------------------------------------------
      SUBROUTINE simsdls ( indic, simext, mct, ydual, funct, grad,
     &                     iwork, dwork)
c-----------------------------------------------------------------------
c
c     INPUT 
c            indic  : = 4, to compute function and gradients     integer
c            simext : entry point of an external subroutine
c                     provided by the user (not used here)
c            mct    : number of variables                        integer
c            ydual  : dual solution                   vector(mct) double
c
c     OUTPUT 
c            funct  : dual function value                         double
c            grad   : gradient value                  vector(mct) double
c
c     WORKSPACE 
c            iwork  : vector (12*nmat+2)                         integer
c                     contents in input :
c                     nmat (size of the problem) in iwork(1)
c                     info (errors of SOLSDLS )  in iwork(2)
c            dwork  : vector                                      double
c                     ( 2*mct + (mct+3)*(nmat*(nmat+1)/2)
c                       + (nmat*(3*nmat+27)) )
c                     contents in input :
c                     vector(mct) b
c                     vectorized symmetric matrix C(nmat*(nmat+1)/2)
c                     mct vectorized symmetric matrices Ai
c                                           mct  * (nmat*(nmat+1)/2)
c     CALL   
c        SOLSDLS : computing the solution of SDLS optimization from the
c                  dual solution
c        NMS     : computing the Frobenius norm of
c                  a vectorized symmetric matrix
c        XV      : computing scalar product of two vectors
c        OPERATA : computing the operator_A : vector(m) result of
c                  scalar product of m Ai constraints matrices with
c                  matrix X ( vectorized version )
c        DV      : computing the difference of 2 vectors
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, mct, info
      DOUBLE PRECISION funct, ydual(*), grad(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables 
      INTEGER nmat, ssm, pb, pc, pai, px, pdw, piw, pinfo, pnmat, pvoa
      DOUBLE PRECISION norm2, yb
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space  : iwork
c     ---------------------------------------- 
      pnmat  = 1
c     pnmat  : pointer for nmat the size of the problem, so 1 more
      pinfo  = pnmat + 1
c     pinfo  : pointer for info (errors of SOLSDLS), so 1 more
      piw    = pinfo + 1
c     piw    : pointer for workspace used by SOLSDLS, so 12*nmat more
c     pinext : pointer for the next iwork array
c              pinext = piw + 12*nmat
c             
c     Total so size of iwork array (12*nmat+2)
c
c     initializations
      info = 0
      nmat = iwork(pnmat)
      ssm = (nmat*(nmat+1)/2)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pb = 1
c     pb    : pointer for the vector B(mct), so mct more
      pc = pb  + mct
c     pc    : pointer for the vectorized symmetric matrix C,
c             so ssm + (nmat*(nmat+1)/2) more
      pai = pc + ssm
c     pai   : pointer for the mct vectorized symmetric matrix Ai
c             (nmat*(nmat+1)/2), so mct*(nmat*(nmat+1)/2) more
      px  = pai +  mct*ssm
c     px    : pointer for the solution x vectorized symmetric matrix
c             so ssm = (nmat*(nmat+1)/2) more
      pvoa= px + ssm
c     pvoa  : pointer for operator_A vector, so mct more
      pdw = pvoa + mct
c     pdw   : pointer for workspace used by SOLSDLS
c             so ( ssm + nmat*(3*nmat+27) ) more
c
c     Total size of dwork array =  mct + ssm + mct*ssm + ssm + mct
c                                + ssm + nmat*(3*nmat+27)
c
c                    = 2*mct + (mct+3)*ssm  + (nmat*(3*nmat+27))
c
c-----------------------------------------------------------------------
c
c     projection x = pK( C+Aiy(i) )
      CALL solsdls ( nmat, dwork(pc), mct, dwork(pai), ydual,
     &               iwork(piw), dwork(pdw),
     &               dwork(px), iwork(pinfo) )
c
c     computing ||x||**2
      CALL NMS ( nmat, dwork(px), norm2 )
      norm2 = norm2*norm2
c
c     computing scalar product ydual by b
      CALL XV ( mct, dwork(pb), ydual, yb )
c
c     computing function value f = -yb + 0.5*||x||**2
      funct = -yb + 0.5*norm2
c
c     computing operator_A = Ai.x
      CALL operata ( mct, nmat, dwork(px), dwork(pai), dwork(pvoa) )
c
c     computing gradient, grad(i) = operator_A - B(i)
      CALL DV ( mct, dwork(pvoa), dwork(pb), grad )
c
      nmat = indic
c
      RETURN
      END
