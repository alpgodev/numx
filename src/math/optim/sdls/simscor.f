c=======================================================================
c
c     SIMSCOR
c
c     Subroutine used by BFGS to compute the value of the dual function
c     and its gradient - specific to obtain a correlation matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE simscor ( indic, simext, mct, ydual, funct, grad,
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
c            iwork  : ( 12*nmat + 1 )                            integer
c                     contents in output :
c                     info (errors of SOLSCOR )  in iwork(1)
c            dwork  : vector ( mct*(5*mct+30) )                   double
c                     contents in input :
c                     vectorized symmetric matrix C(nmat*(nmat+1)/2)
c     CALL   
c        SOLSCOR : computing the solution of SDLS optimization from the
c                  dual solution
c                  specific to obtain a correlation matrix
c        NMS     : computing the Frobenius norm of
c                  a vectorized symmetric matrix
c        CMDSV   : converting a vectorized symmetric diagonal matrix
c                  in vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, mct
      DOUBLE PRECISION funct, ydual(*), grad(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER nmat, ssm, i, pc, px, pdw, piw, pinfo, pvoa
      DOUBLE PRECISION norm2, yb
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pinfo  = 1
c     pinfo  : pointer for info (errors of SOLSCOR), so 1 more
      piw    = pinfo + 1
c     piw    : pointer for workspace used by SOLSCOR, so 12*nmat more
c     pinext : pointer for the next iwork array
c              pinext = piw + 12*nmat
c
c     Total size of iwork array (12*nmat+1)
c
c     initializations
      nmat = mct
      ssm = (nmat*(nmat+1)/2)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pc = 1
c     pc    : pointer for the vectorized symmetric matrix C,
c             so ssm = (nmat*(nmat+1)/2) more
      px  = pc +  ssm
c     px    : pointer for the solution x vectorized symmetric matrix
c             so ssm = (nmat*(nmat+1)/2) more
      pvoa= px + ssm
c     pvoa  : pointer for operator_A vector, so mct more
      pdw = pvoa + mct
c     pdw   : pointer for workspace used by SOLSCOR
c             so ( ssm + nmat*(3*nmat+27) ) more
c     pnext :  pointer for the next work array
c              pnext = pdw + ( ssm + nmat*(3*nmat+27) )
c
c              dwtot = ssm + ssm + ssm + mct
c                    +  ( ssm + nmat*(3*nmat+27) )
c                    = mct + 4*ssm + (nmat*(3*nmat+27))
c              for correlation mct = nmat, so
c                    = mct + 4*ssm + (mct*(3*mct+27))
c                    = 4*ssm + (mct*(4*mct+28))
c                    = 4*(mct*(mct+1)/2) + (mct*(3*mct+28))
c                    = mct * ( (4*(mct+1)/2) + (3*mct+28) )
c                    = mct * ( (2*(mct+1) + (3*mct+28) )
c                    = mct * ( (2*mct+2) + (3*mct+28) )
c     Total size of dwork array = mct * (5*mct+30)
c
c-----------------------------------------------------------------------
c
c     projection x = pK( C+Aiy(i) )
      CALL solscor ( nmat, dwork(pc), ydual,
     &               iwork(piw), dwork(pdw),
     &               dwork(px), iwork(pinfo) )
c
c     computing ||x||**2
      CALL NMS ( nmat, dwork(px), norm2 )
      norm2 = norm2*norm2
c
c     computing scalar product ydual by b (vector of 1 for correlation)
      yb = 0.
      DO i = 1,mct
         yb = yb + ydual(i)
      ENDDO
c
c     computing function value 
c     f = -y'b + 0.5||pK(C+Ay)||**2
c       = -yb + 0.5*||x**2||
      funct = -yb + 0.5*norm2
c
c     computing operator_A = Ai.x
c     (vector of diagonal(X) for correlation)
c     Converting the diagonal of a symmetric matrix(n*(n+1)/2) in vector(n)
      CALL CMDSV ( nmat, dwork(px), dwork(pvoa) )
c
c     computing gradient, grad(i) = operator_A - b(i)
      DO i = 1,mct
         grad(i) = dwork(pvoa+i-1) - 1
      ENDDO
c
      nmat = indic
      RETURN
      END
