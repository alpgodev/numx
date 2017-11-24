c=======================================================================
c
c     BFGSBOXS                                   
c
c     Quasi-Newton optimizer with BFGS method.
c
c     min f(x) subject to: binf <= x <= bsup
c
c-----------------------------------------------------------------------
      SUBROUTINE bfgsboxs ( simul, simext, nvar, xvar, epsbfg,
     &                      binf, bsup, iwork, dwork,
     &                      xopt, funct, grad, totitr, totsim, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            SIMUL  : entry point of objective function and gradient
c            SIMEXT : entry point of external subroutine
c            nvar   : number of variables                        integer
c            xvar   : initial value of the variables (nvar)       double
c            epsbfg : precision stop test                         double
c            binf   : inferior bounds (nvar)                      double
c            bsup   : superior bounds (nvar)                      double
c
c
c     WORKSPACE 
c            iwork  : 2*nvar+1 + space for SIMUL                 integer
c            dwork  : nvar*(nvar+11)/2 + space for SIMUL          double 
c                    
c     OUTPUT 
c            xopt   : optimal value(s) at the minimum (nvar)      double
c            funct  : objective function at the minimum           double
c            grad   : gradients value(s) at the minimum (nvar)    double
c            totitr : number of iterations                       integer
c            totsim : number of calls of SIMUL                   integer
c            info   : diagnostic argument                        integer
c
c     CALL   
c            SIMUL   : subroutine who computes the function and gradient
c            QNBFGS  : quasi-Newton optimizer with BFGS method (optim)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simul, simext
c
c     i/o arguments
      INTEGER nvar, totitr, totsim, info
      DOUBLE PRECISION epsbfg, funct
      DOUBLE PRECISION xvar(*), binf(*), bsup(*), grad(*), xopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER imp, ion2qn, indic, minmax, mode, i, pdxmin, piwb, pdwb,
     &        piws, pdws
      REAL*4 rzs(1)
      DOUBLE PRECISION df1, epsabs, epsdxm, normg, normin
      PARAMETER ( minmax=300, normin=1.d-6, epsdxm=1.d-10 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info   = 0
      imp    = 0
      ion2qn = 6
      indic  = 4
      info   = 0
      mode   = 1
      totitr =  max(minmax, 2*nvar)
      totsim = totitr*3
c
c     pointers for integer workspace : iwork
c     --------------------------------------
      piwb = 1
c     piwb : pointer for n2qn1 who needs (2*nvar + 1)
c             so (2*nvar + 1) more
      piws = piwb + ( 2*nvar + 1 )
c     piws : pointer for SIMUL
c
c     pointers for double precision  workspace : dwork
c     ------------------------------------------------
      pdxmin = 1
c     pdxmin : pointer for dxmin, so (nvar) more
      pdwb = pdxmin + nvar
c     pdwb : pointer for n2qn1 who needs (nvar*(nvar+9)/2)
c             so (nvar*(nvar+9)/2) more
      pdws = pdwb +  nvar*(nvar+9)/2
c     pdws : pointer for SIMUL
c
c-----------------------------------------------------------------------
c
c     first call to initialize function and gradient values
      CALL simul ( indic, simext, nvar, xvar, funct, grad,
     &             iwork(piws), dwork(pdws) )
c
c     computing the precision stop criter
c     ( when mean of projecting gradient < epsabs )
      epsabs = epsbfg * nvar
c      epsabs = epsbfg * sqrt(float(nvar))
c      epsabs = epsbfg
c
c     computing df1 : estimation of the criter decrease
c                     at the first iteration
c                     must be nearly abs( f(xvar)-f(xopt) ) )
      df1 = abs(funct) / 10.
c
c     computing the precision on variables stop criter
      DO i=1,nvar
         normg = max(abs(grad(i)), normin)
         dwork(pdxmin+i-1) = ( df1* epsdxm ) / normg
         xopt(i) = xvar(i)
c         dwork(pdxmin+i-1) = ( df1* epsbfg ) / normg
c         dwork(pdxmin+i-1) = ( df1* epsabs ) / normg
c         dwork(pdxmin+i-1) = sqrt (( df1* epsdxm ) / normg)
      ENDDO
c
c     non linear optimization - BFGS Solver
       CALL qnbfgs ( simul, simext, nvar, xopt,
     &             funct, grad, dwork(pdxmin), df1, epsabs,
     &             imp, ion2qn ,mode, totitr, totsim, binf, bsup,
     &             iwork(piwb), dwork(pdwb),
     &             iwork(piws), rzs, dwork(pdws) )
c
c     error codes
      IF (mode.le.0) info = -1101
      IF (mode.eq.1) info = 0
      IF (mode.eq.2) info = -1102
      IF (mode.eq.3) info = -1103
      IF (mode.eq.4) info =  1104
      IF (mode.eq.5) info =  1105    
      IF (mode.eq.6) info =  1106
      IF (mode.eq.7) info = -1107
c
      RETURN
      END
