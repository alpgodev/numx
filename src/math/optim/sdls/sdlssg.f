c=======================================================================
c
c     SDLSSG                                
c
c     Semi-Definite Least Square optimization
c     general version with equal and inequal constraints
c     ( the constraints matrices Ai and Bj are symmetric )
c
c-----------------------------------------------------------------------
      SUBROUTINE sdlssg ( n, cmat, mcte, acte, bcte,
     &                    mcti, acti, lcti, ucti, epsbfg,
     &                    iwork, dwork, xmat, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix (nmat*nmat) to optimize
c                     vector(nmat*(nmat+1)/2)                     double
c            mcte   : number of constraints                      integer
c            acte   : mcte Ai symmetric matrices(nmat*nmat)
c                     of constraints
c                     vector( mcte*(nmat*(nmat+1)/2) )            double
c            bcte   : vector(mcte) of constraints                 double
c            mcti   : number of inequal constraints              integer
c            acti   : mcti Bj symmetric matrices(nmat*nmat) of
c                     inequal constraints  vector(mcti*nmat*nmat)
c                     vector( mcti*(nmat*(nmat+1)/2) )            double
c            lcti   : lower bounds of inequal constraints
c                     vector(mcti)                                double
c            ucti   : upper bounds of inequal constraints
c                     vector(mcti)                                double
c            epsbfg : BFGS precision stop test                    double
c
c     WORKSPACE 
c            iwork  : ( 2*mcte + 4*mcti + 12*nmat + 5 )          integer
c            dwork  : vector                                      double
c                     (  (mcte+2*mcti)*(mcte+2*mcti+23)/2)
c                        + (mcte+mcti+3)*(nmat*(nmat+1)/2)
c                        + (nmat*(3*nmat+27))  )
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
c        SIMSDLSG   : subroutine used by BFGS to compute the value of the
c                     dual function and its gradient
c        SOLSDLSG   : computing the solution of SDLSG optimization from the
c                     dual solution   x = pK( C + Ai.y(i) ) 
c        PROJECTS   : projection of a symmetric matrix
c                     on the cone of semidefinite matrices
c                     version with vector (n*(n+1)/2) in input/output
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simsdlsg, gesterr
c
c     i/o arguments
      INTEGER n, mcte, mcti, info
      DOUBLE PRECISION epsbfg, cmat(*), xmat(*), acte(*), bcte(*), 
     &                 acti(*), lcti(*), ucti(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mdual, ssm, totitr, totsim, piwbf, pnmat, pmcte, pmcti, 
     &        pinfo, pwbf, pc, psa, pb, psb, pl, pu, pydual, pbinf, 
     &        pbsup, pgrad, pisim, psim, pisol, psol, pipro, ppro, 
     &        pbilam, pbimu, pbigam, pbslam, pbsmu, pbsgam, infotmp
      DOUBLE PRECISION infini, dzero, funct
      PARAMETER ( infini = 1.E20, dzero = 0.0E0 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info  = 0
      infotmp = 0
      ssm   = n*(n+1)/2
      mdual = mcte + 2*mcti
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piwbf = 1
c     piwbf : pointer for BFGSBOX who needs (2*mdual + 1)
c             so (2*mdual + 1) more
      pisim  = piwbf + (2*mdual + 1)
c     pisim : pointer for the simulator SIMSDLSG who needs (12*n+4)
      pnmat = pisim
c     pnmat : pointer for value nmat the size of the problem, so 1 more
      pmcte  = pnmat +1
c     pmcte  : pointer for number of equal constraints mcte, so 1 more
      pmcti  = pmcte +1
c     pmcti  : pointer for number of equal constraints mcti, so 1 more
      pinfo  = pmcti + 1
c     pinfo : pointer for info (errors of simul), so 1 more
      pisol= 1
c     pisol : pointer for SOLSDLSG integer workspace
c             SOLSDLSG needs always less place than SIMSDLSG
c             and used all now free space
      pipro= 1
c     pipro : pointer for PROJECTS integer workspace
c             PROJECTS also needs always less place than SIMSDLSG
c
c     Total size of dwork array : (2*mdual+1) + (12*n+4)
c                               = (2*mdual+12*n+5)
c                               = (2*(mcte+2*mcti)+12*n+5)
c                               = (2*mcte+4*mcti+12*n+5)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pydual = 1
c     pydual: pointer for initial value of ydual, so mdual more
      pbinf = pydual + mdual
c     pbinf : pointer for inferior bounds binf, so mdual more
      pbsup = pbinf + mdual
c     pbsup : pointer for superior bounds bsup, so mdual more
      pgrad = pbsup + mdual
c     pgrad : pointer for gradient at solution, so mdual more
c
      pwbf = pgrad + mdual
c     pwbf  : pointer for BFGSSDLS who needs ( mdual*(mdual+11)/2 )
c             so ( mdual*(mdual+11)/2 ) more
c
      psim = pwbf + mdual*(mdual+11)/2
c     psim  : pointer for the simulator SIMSDLSG who needs
c                     ( (3+mcte+mcti)*(n*(n+1)/2)
c                            + 2*mdual + n*(3*n+27) )
      pc = psim
c     pc    : pointer for the vectorized symmetric matrix C,
c             communication with SIMSDLSG,
c             so ssm =  (nmat*(nmat+1)/2) more
      psa = pc + ssm
c     psa   : pointer for the mcte vectorized symmetric matrix
c             Ai(n*(n+1)/2), communication with SIMSDLSG,
c             so mcte*(n*(n+1)/2) more
      pb = psa + mcte*ssm
c     pb    : pointer for the vector b(mcte), communication
c             with SIMSDLSG, so mcte more
      psb = pb + mcte
c     psb   : pointer for the mcti vectorized symmetric matrix
c             Bj(n*(n+1)/2), communication with SIMSDLSG,
c             so mcti*(n*(n+1)/2) more
      pl = psb + mcti*ssm
c     pl    : pointer for the vector low(mcti), communication
c             with SIMSDLSG, so mcti more
      pu = pl + mcti
c     pu    : pointer for the vector low(mcti), communication
c             with SIMSDLSG, so mcti more
      psol= pwbf
c     psol  : pointer for SOLSDLSG workspace
c             SOLSDLSG needs always less place than SIMSDLSG because
c             SIMSDLSG call also SOLSDLSG
c             and used now free space of BFGSSDLS + SIMSDLSG
      ppro= 1
c     ppro  : pointer for PROJECTS workspace 
c             PROJECTS also needs always less place than SIMSDLSG
c
c     Totam size of dwork array : 4*mdual + mdual*(mdual+11)/2
c                               + (3+mcte+mcti)*(n*(n+1)/2)
c                               + 2*mdual + n*(3*n+27)
c                         = 6*(mcte+2*mcti)
c                         + (mcte+mcti)*(mcte+2*mcti+11)/2
c                         + (mcte+mcti+3)*(n*(n+1)/2)
c                         + (n*(3*n+27))
c
c                         = (mcte+2*mcti)*(mcte+2*mcti+23)/2)
c                         + (mcte+mcti+3)*(n*(n+1)/2)
c                         + (n*(3*n+27))
c
c     pointers for bounds : at dwork(pbinf) and dwork(pbsup)
c     ------------------------------------------------------
      pbilam = pbinf
      pbimu  = pbilam + mcte
      pbigam = pbimu + mcti
      pbslam = pbsup
      pbsmu  = pbslam + mcte
      pbsgam = pbsmu + mcti
c
c-----------------------------------------------------------------------
c     case with constraints
c
      IF (mdual .NE. 0) THEN
c
c     initialization of ydual
      CALL IVX ( mdual, dwork(pydual), dzero )
c
c     initializations of bounds
      CALL IVX ( mcte, dwork(pbilam), -infini )
      CALL IVX ( mcte, dwork(pbslam), infini )
      CALL IVX ( mcti, dwork(pbimu), dzero )
      CALL IVX ( mcti, dwork(pbsmu), infini )
      CALL IVX ( mcti, dwork(pbigam), dzero )
      CALL IVX ( mcti, dwork(pbsgam), infini )
c
c     construction of work vectors (communication with simsdlsg)
      iwork(pnmat) = n
      iwork(pmcte) = mcte
      iwork(pmcti) = mcti
      iwork(pinfo) = 0
c     matrix C (n*(n+1)/2)
      CALL YV ( ssm, cmat, dwork(pc) )
c     matrices of equal constraints  mcte*(n*(n+1)/2)
      CALL YV ( mcte*ssm, acte, dwork(psa) )
c     vector of equal constraints b (mcte)
      CALL YV ( mcte, bcte, dwork(pb) )
c     matrices of inequal constraints  mcti*(n*(n+1)/2)
      CALL YV ( mcti*ssm, acti, dwork(psb) )
c     vector of low inequal constraints (mcti)
      CALL YV ( mcti, lcti, dwork(pl) )
c     vector of up inequal constraints (mcti)
      CALL YV ( mcti, ucti, dwork(pu) )
c
c     optimization BFGS SDLS 
c
c     (the second parameter is not used here, we can replace GESTERR
c      by an other name of subroutine)
      CALL bfgssdls ( simsdlsg, gesterr, mdual, dwork(pydual), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pwbf),
     &                funct, dwork(pgrad), totitr, totsim, infotmp)
c
c     construction of the matrix solution from dual solution
      CALL solsdlsg ( n, cmat, mcte, acte,
     &                mcti, acti, dwork(pydual),
     &                iwork(pisol), dwork(psol),
     &                xmat, info )
      IF (info .LT. 0) RETURN
c
c     case with no constraints
c     projection on sdp matrix cone
      ELSE
         CALL PROJECTS ( n, cmat, iwork(pipro),dwork(ppro), xmat, info )
      ENDIF
      info = infotmp
c
      RETURN
      END
