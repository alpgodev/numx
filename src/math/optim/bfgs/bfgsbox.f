c=======================================================================
c
c     BFGSBOX                                    
c
c     Quasi-Newton optimizer with BFGS method - Expert version
c
c     min f(x) s.t. binf <= x <= bsup  
c
c-----------------------------------------------------------------------
      SUBROUTINE bfgsbox ( SIMUL, SIMEXT, nvar, xvar, dxmin, df1,
     &                     epsabs, mode, maxitr, maxsim, binf, bsup,
     &                     iwork, dwork,
     &                     xopt, funct, grad, totitr, totsim, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            SIMUL  : entry point of objective function and gradient
c            SIMEXT : entry point of external subroutine
c            nvar   : number of variables                        integer
c            xvar   : initial value of the variables (nvar)       double
c            dxmin  : precision on the variables (nvar)           double
c            df1    : estimation of the first decrease of the criter
c                     used if mode=1 computed if =<0 
c                     must be nearly abs( f(xvar)-f(xopt) ) )     double
c            epsabs : precision stop test                         double
c            mode   : type of initialization                     integer
c                     = 1 : nothing to initialize
c            maxitr : maximum number of iteration                integer
c            maxsim : maximum number of calls of SIMUL           integer
c            binf   : lower bounds (nvar)                         double
c            bsup   : upper bounds (nvar)                         double
c
c     WORKSPACE 
c            iwork  : ( 2*nvar+1 ) + space for SIMUL             integer
c            dwork  : ( nvar*(nvar+9)/2) ) + space for SIMUL      double
c
c     OUTPUT 
c            xopt   : optimal points at the minimum (nvar)        double
c            funct  : objective function at the minimum           double
c            grad   : gradient(s) at the minimum (nvar)           double
c            totitr : number of iteration(s)                     integer
c            totsim : number of call(s) of SIMUL                 integer
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
      INTEGER nvar, mode, maxitr, maxsim, totitr, totsim, info
      DOUBLE PRECISION funct, df1, epsabs, xopt(*)
      DOUBLE PRECISION xvar(*), dxmin(*), binf(*), bsup(*), grad(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER imp, ion2qn, indic, piwb, pdwb, piws, pdws, i
      DOUBLE PRECISION epsn2q
      REAL*4 rzs(1)
      INTEGER iosort
      PARAMETER ( iosort=66 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info   = 0
      imp    = 0
      ion2qn = 6
      indic  = 4
      info   = 0
      totitr = maxitr
      totsim = maxsim
      epsn2q = epsabs
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
      pdwb = 1
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
c     computing df1 if not provided
      IF ( df1.le.0) df1 = abs(funct) / 10.
      
c
c     computing the precision on variables stop criter
      DO i=1,nvar
         xopt(i) = xvar(i)
      ENDDO
      
c
c     non linear optimization - BFGS
      CALL qnbfgs ( simul, simext, nvar, xopt,
     &             funct, grad, dxmin, df1, epsn2q,
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
