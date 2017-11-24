c=======================================================================
c
c     SIMSDLSG
c
c     Subroutine used by BFGS to compute the value of the dual function
c     and its gradient
c
c-----------------------------------------------------------------------
      SUBROUTINE simsdlsg ( indic, simext, mdual, ydual, funct, grad,
     &                      iwork, dwork)
c-----------------------------------------------------------------------
c
c     INPUT 
c            indic  : = 4, to compute function and gradients     integer
c            simext : entry point of an external subroutine
c                     provided by the user (not used here)
c            mdual  : number of variables                        integer
c            ydual  : dual solution                 vector(mdual) double
c
c     OUTPUT 
c            funct  : dual function value                         double
c            grad   : gradient value                vector(mdual) double
c
c     WORKSPACE 
c            iwork  : ( 12*nmat + 4 )                            integer
c                     contents in input :
c                     nmat (size of the problem) in iwork(1)
c                     mcte (nb of equal constraints) in iwork(2)
c                     mcti (nb of inequal constraints) in iwork(3)
c                     info (errors of SOLSDLS )  in iwork(4)
c            dwork  : vector                                      double
c                     ( (3+mcte+mcti)*(nmat*(nmat+1)/2)
c                            + 2*mdual + nmat*(3*nmat+27) )
c                     contents in input :
c                     vectorized symmetric matrix C(nmat*(nmat+1)/2)
c                     mcte vectorized symmetric matrices Ai of equal
c                              constraints    mcte * (nmat*(nmat+1)/2)
c                     vector(mcte) b of equal constraints
c                     mcti vectorized symmetric matrices Bj of inequal
c                              constraints    mcti * (nmat*(nmat+1)/2)
c                     vector(mcti) low bounds of inequal constraints
c                     vector(mcti) up bounds of inequal constraints
c
c     CALL   
c        SOLSDLSG : computing the solution of SDLSG optimization
c                   from the dual solution
c        NMS     : computing the Frobenius norm of
c                  a vectorized symmetric matrix
c        XV       : computing scalar product of two vectors
c        OPERATA  : computing the operator_A : vector(m) result of
c                   scalar product of m Ai constraints matrices with
c                   matrix X ( vectorized version )
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c        DV       : computing the difference of 2 vectors
c        SV      : computing the sum of 2 vectors
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, mdual
      DOUBLE PRECISION funct, ydual(*), grad(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER nmat, mcte, mcti, ssm, pnmat, pmcte, pmcti, pinfo, piw, 
     &        pc, psa, pb, psb, pl, pu, px, pdw, pvoc, plam, pgam, pmu,
     &        pgb, pgl, pgu
      DOUBLE PRECISION norm2, yb, yl, yu, yres, dmun
      PARAMETER ( dmun = -1.0E0 )
c
c-----------------------------------------------------------------------
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pnmat  = 1
c     pnmat  : pointer for nmat the size of the problem, so 1 more
      pmcte  = pnmat +1
c     pmcte  : pointer for number of equal constraints mcte, so 1 more
      pmcti  = pmcte +1
c     pmcti  : pointer for number of equal constraints mcti, so 1 more
      pinfo  = pmcti + 1
c     pinfo  : pointer for info (errors of SOLSDLS), so 1 more
      piw    = pinfo + 1
c     piw    : pointer for workspace used by SOLSDLS, so 12*nmat more
c     pinext : pointer for the next iwork array
c              pinext = piw + 12*nmat
c              so size of iwork array (12*nmat+4)
c
c     initializations
      nmat = iwork(pnmat)
      mcte = iwork(pmcte)
      mcti = iwork(pmcti)
      ssm = (nmat*(nmat+1)/2)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pc = 1
c     pc    : pointer for the vectorized symmetric matrix C,
c             so ssm = (nmat*(nmat+1)/2) more
      psa = pc + ssm
c     psa   : pointer for the mcte vectorized symmetric matrix Ai
c             (nmat*(nmat+1)/2), so mcte*(nmat*(nmat+1)/2) more
      pb = psa + mcte*ssm
c     pb    : pointer for the vector B(mcte), so mcte more
      psb = pb + mcte
c     psb   : pointer for the mcti vectorized symmetric matrix
c             Bj(nmat*(nmat+1)/2), so mcti*(nmat*(nmat+1)/2) more
      pl = psb + mcti*ssm
c     pl    : pointer for the vector low(mcti), communication
c             with SIMSDLSG, so mcti more
      pu = pl + mcti
c     pu    : pointer for the vector up(mcti), communication
c             with SIMSDLSG, so mcti more
      px = pu + mcti
c     px    : pointer for the solution x vectorized symmetric matrix
c             so ssm = (nmat*(nmat+1)/2) more
      pvoc= px + ssm
c     pvoc  : pointer for operator_C vector, so mdual more
      pdw = pvoc + mdual
c     pdw   : pointer for workspace used by SOLSDLS
c             so ( ssm + nmat*(3*nmat+27) ) more
c
c     Total size of work array = ssm + mcte*ssm + mcte
c                              + mcti*ssm + mcti + mcti
c                              + ssm + mdual
c                              + ssm + nmat*(3*nmat+27)
c
c                        = (3+mcte+mcti)*ssm + mcte + 2*mcti + mdual
c                                + nmat*(3*nmat+27)
c                          with mdual = mcte + 2*mcti
c
c                        = (3+mcte+mcti)*(nmat*(nmat+1)/2)
c                            + 2*mdual + nmat*(3*nmat+27)
c
c     pointers for ydual vector(mdual)
c     --------------------------------
      pgb = 1
      pgl = pgb + mcte
      pgu = pgl + mcti
c
c     pointers for operator_C vector(mdual) : at dwork(pvoc)
      plam = pvoc
      pmu = plam + mcte
      pgam = pmu + mcti
c
c-----------------------------------------------------------------------
c
c     projection x = pK( C+Aiy(i) )
      CALL solsdlsg ( nmat, dwork(pc), mcte, dwork(psa),
     &                mcti, dwork(psb), ydual,
     &                iwork(piw), dwork(pdw),
     &                dwork(px), iwork(pinfo) )
c
c     computing ||x||**2
      CALL NMS ( nmat, dwork(px), norm2 )
      norm2 = norm2*norm2
c
c     computing scalar product ydual by constraint values vector
      CALL XV ( mcte, dwork(pb), ydual(pgb), yb )
      CALL XV ( mcti, dwork(pl), ydual(pgl), yl )
      CALL XV ( mcti, dwork(pu), ydual(pgu), yu )
      yres = yb + yl - yu
c
c     computing function value f = -yres + 0.5*||x**2||
      funct = -yres + 0.5*norm2
c
c     computing operator_C (vector(mdual)) = Ci.x = Ai.x;Bj.x;-Bj.x
      CALL operata ( mcte, nmat, dwork(px), dwork(psa), dwork(plam) )
      CALL operata ( mcti, nmat, dwork(px), dwork(psb), dwork(pmu) )
      CALL PVX ( mcti, dwork(pmu), dmun, dwork(pgam) )
c
c     computing gradient, grad(k) = operator_C - (b(i);low(j);-up(j))
      CALL DV ( mcte, dwork(plam), dwork(pb), grad(pgb) )
      CALL DV ( mcti, dwork(pmu), dwork(pl), grad(pgl) )
      CALL SV ( mcti, dwork(pgam), dwork(pu), grad(pgu) )
c
      nmat = indic
      RETURN
      END
