c=======================================================================
c
c     SDLSS                                
c
c     Semi-Definite Least Square optimization
c     ( the constraints matrices Ai are symmetric )
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlss ( nmat, cmat, mct, aict, bct, epsbfg,
     &                   iwork, dwork, xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix (nmat*nmat) to optimize
c                     vector(nmat*(nmat+1)/2)                     double
c            mct    : number of constraints                      integer
c            aict   : mct Ai symmetric matrices(nmat*nmat)
c                     of constraints
c                     vector( mct*(nmat*(nmat+1)/2) )             double
c            bct    : vector(mct) of constraints                  double
c            epsbfg : BFGS precision stop test                    double
c
c     WORKSPACE 
c            iwork  : 2*mct + 12*nmat + 3                        integer
c            dwork  : ( mct*(mct+23)/2)
c                    + (mct+3)*(nmat*(nmat+1)/2)
c                    + (nmat*(3*nmat+27)) )                       double 
c
c     OUTPUT 
c            xmat   : matrix solution (nmat*(nmat+1)/2)           double
c            info   : diagnostic argument                        integer
c
c     CALL   
c        IVX        : initialization at a scalar of a vector
c        YV         : copy a vector in a vector
c        BFGSSDLS   : quasi-Newton optimizer with BFGS method
c                     ( short call version for SDLS )
c        SIMSDLS    : subroutine used by BFGS to compute the value of the
c                     dual function and its gradient
c        SOLSDLS    : computing the solution of SDLS optimization from the
c                     dual solution   x = pK( C + Ai.y(i) ) 
c        PROJECTS   : projection of a symmetric matrix
c                     on the cone of semidefinite matrices
c                     version with vector (n*(n+1)/2) in input/output
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simsdls, gesterr
c
c     i/o arguments
      INTEGER nmat, mct, info
      DOUBLE PRECISION epsbfg, cmat(*), aict(*), bct(*), xmat(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER ssm, totitr, totsim, piwbf, pnmat, pinfo, pwbf, pb, pc, 
     &        pai, pydual, pbinf, pbsup, pgrad, pisim, psim, pisol, 
     &        psol, pipro, ppro, infotmp
      DOUBLE PRECISION infini, dzero, funct
      PARAMETER ( infini = 1.E20, dzero = 0.E0 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      infotmp = 0
      ssm  = nmat*(nmat+1)/2
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piwbf = 1
c     piwbf : pointer for BFGSBOX who needs (2*mct + 1)
c             so (2*mct + 1) more
      pisim  = piwbf + (2*mct + 1)
c     pisim : pointer for the simulator SIMSDLS who needs (12*nmat+2)
      pnmat = pisim
c     pnmat : pointer for value nmat the size of the problem, so 1 more
      pinfo = pnmat + 1
c     pinfo : pointer for info (errors of simul), so 1 more
      pisol= 1
c     pisol : pointer for SOLSDLS integer workspace
c             SOLSDLS needs always less place than SIMSDLS
c             and used all now free space
      pipro= 1
c     pipro : pointer for PROJECTS integer workspace
c             SOLSDLS also needs always less place than SIMSDLS
c
c     Total size of dwork array : (2*mct+1) + (12*nmat+2)
c                               = (2*mct+12*nmat+3)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pydual = 1
c     pydual: pointer for initial value of ydual, so mct more
      pbinf = pydual + mct
c     pbinf : pointer for inferior bounds binf, so mct more
      pbsup = pbinf + mct
c     pbsup : pointer for superior bounds bsup, so mct more
      pgrad = pbsup + mct
c     pgrad : pointer for gradient at solution, so mct more
c
      pwbf = pgrad + mct
c     pwbf  : pointer for BFGSSDLS who needs ( mct*(mct+11)/2 )
c             so ( mct*(mct+11)/2 ) more
c
      psim = pwbf +  mct*(mct+11)/2
c     psim  : pointer for the simulator SIMSDLS who needs
c                     ( (3+mct)*(nmat*(nmat+1)/2)
c                            + 2*mct + nmat*(3*nmat+27) )
      pb = psim
c     pb    : pointer for the vector B(mct), communication with SIMSDLS,
c             so mct more
      pc = pb  + mct
c     pc    : pointer for the vectorized symmetric matrix C,
c             communication with SIMSDLS,
c             so ssm =  (nmat*(nmat+1)/2) more
      pai = pc + ssm
c     pai   : pointer for the mct vectorized symmetric matrix
c             Ai(nmat*(nmat+1)/2), communication with SIMSDLS,
c             so mct*(nmat*(nmat+1)/2) more
      psol= pwbf
c     psol  : pointer for SOLSDLS workspace
c             SOLSDLS needs always less place than SIMSDLS because
c             SIMSDLS call also SOLSDLS
c             and used now free space of BFGSSDLS + SIMSDLS
      ppro= 1
c     ppro  : pointer for PROJECTS workspace 
c             PROJECTS also needs always less place than SIMSDLS
c
c     Total size of dwork array : 4*mct + mct*(mct+11)/2
c                         + 2*mct + (mct+3)*(nmat*(nmat+1)/2)
c                         + nmat*(3*nmat+27)
c     
c                         = mct*(mct+11)/2) + 6*mct
c                         + (mct+3)*(nmat*(nmat+1)/2)
c                         + (nmat*(3*nmat+27))
c
c                         = mct*(mct+23)/2)
c                         + (mct+3)*(nmat*(nmat+1)/2)
c                         + (nmat*(3*nmat+27))
c     
c-----------------------------------------------------------------------
c
      IF (mct .NE. 0) THEN
c
c     initialization of ydual
      CALL IVX ( mct, dwork(pydual), dzero )
c
c     initializations of bounds
      CALL IVX ( mct, dwork(pbinf), -infini )
      CALL IVX ( mct, dwork(pbsup), infini )
c
c     construction of work vectors (communication with simsdls)
      iwork(pnmat) = nmat
      iwork(pinfo) = 0
c     matrix C (nmat*(nmat+1)/2)
      CALL YV ( ssm, cmat, dwork(pc) )
c     matrices of equal constraints Ai mct*(nmat*(nmat+1)/2)
      CALL YV ( mct*ssm, aict, dwork(pai) )
c     vector of equal constraints b (mct)
      CALL YV ( mct, bct, dwork(pb) )
c
c     optimization BFGS-SDLS
c
c     (the second parameter is not used here, we can replace GESTERR
c      by an other name of subroutine)
      CALL bfgssdls ( simsdls, gesterr, mct, dwork(pydual), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pwbf),
     &                funct, dwork(pgrad), totitr, totsim, infotmp)
c
c     construction of the matrix solution from dual solution
      CALL solsdls ( nmat, cmat, mct, aict, dwork(pydual),
     &               iwork(pisol), dwork(psol),
     &               xmat, info )
      IF (info .LT. 0) RETURN
c
c     case with no constraints : we do only a projection
c                                on the cone of semidefinite matrices
      ELSE
         CALL PROJECTS ( nmat, cmat, iwork(pipro), dwork(ppro),
     &                   xmat, info )
         IF (info .LT. 0) RETURN
      ENDIF
      
      info = infotmp
      
      RETURN
      END
