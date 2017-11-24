c=======================================================================
c
c     subroutine  BFGSMVOL
c
c     Quasi-Newton optimizer with BFGS method
c     (short call version for MULTIVOL)
c
c-----------------------------------------------------------------------
      SUBROUTINE bfgsmvol ( simul, simext, nvar, xvar, epsbfg,
     &                      binf, bsup, iwork, dwork,
     &                      funct, grad, totitr, totsim, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       simul  : entry point of the subroutine who computes
c                the function and the gradient
c       simext : entry point of an external subroutine
c       nvar   : number of variables                            integer
c       xvar   : initial value of the variables vector(nvar)     double
c       epsbfg : precision stop test                             double
c       binf   : inferior bounds on variables   vector(nvar)     double
c       bsup   : superior bounds on variables   vector(nvar)     double
c
c
c     WORKSPACE 
c       iwork  : 2*nvar+1                                        integer
c                + space for SIMUL
c       dwork  : nvar*(nvar+11)/2 
c                + space for SIMUL                                double
c                    
c     OUTPUT 
c       xvar   : variables value at the minimum,vector(nvar)      double
c       funct  : minimum of the function                          double
c       grad   : gradients value at the minimum,vector(nvar)      double
c       totitr : total of iterations                             integer
c       totsim : total of calls of SIMUL                         integer
c       info   : = 0 successful exit                             integer
c
c     CALL   
c        SIMUL   : subroutine who computes the function and gradient
c        QNBFGS  : quasi-Newton optimizer with BFGS method (optim)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simul, simext
c
c     arguments i/o
      INTEGER nvar, totitr, totsim, info
      DOUBLE PRECISION epsbfg, funct
      DOUBLE PRECISION xvar(*), binf(*), bsup(*), grad(*)
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
c
      integer iosort
      parameter ( iosort=6 )
c
      PARAMETER ( minmax=300, normin=1.d-6, epsdxm=1.d-10 )
c     
c     intrinsic functions
      INTRINSIC MAX, ABS
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
      totitr =  MAX(minmax,2*nvar)
      totsim = totitr*3
c
c     pointers for integer workspace : vector iwork
c     ---------------------------------------------
      piwb = 1
c     piwb : pointer for QNBFGS who needs (2*nvar + 1)
c             so (2*nvar + 1) more
      piws = piwb + ( 2*nvar + 1 )
c     piws : pointer for SIMUL
c
c     pointers for double precision  workspace : vector dwork
c     -------------------------------------------------------
      pdxmin = 1
c     pdxmin : pointer for dxmin, so (nvar) more
      pdwb = pdxmin + nvar
c     pdwb : pointer for QNBFGS who needs (nvar*(nvar+9)/2)
c             so (nvar*(nvar+9)/2) more
      pdws = pdwb +  (nvar*(nvar+9))/2
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
c
c     computing df1 : estimation of the criter decrease
c                     at the first iteration
c                     must be nearly abs( f(xvar)-f(xopt) ) )
      df1 = ABS(funct) / 10.
c
c     computing the precision on variables stop criter

      normg = normin
      DO i = 1,nvar
         normg = MAX(ABS(grad(i)),normg)
      ENDDO
      DO i = 1,nvar
         dwork(pdxmin+i-1) = ( df1*epsdxm ) / normg
      ENDDO
c
c     optimization Quasi-Newton
      CALL qnbfgs ( simul, simext, nvar, xvar,
     &             funct, grad, dwork(pdxmin), df1, epsabs,
     &             imp, ion2qn ,mode, totitr, totsim, binf, bsup,
     &             iwork(piwb), dwork(pdwb),
     &             iwork(piws), rzs, dwork(pdws) )
c
c     error management
      IF (mode.le.0) info = -1101   ! error in simulator function
      IF (mode.eq.1) info = 0       ! all is right
      IF (mode.eq.2) info = -1102   ! an input parameter is badly init.
      IF (mode.eq.3) info = -1103   ! matrix no sdp
      IF (mode.eq.4) info =  1104   ! max. number of iterations
      IF (mode.eq.5) info =  1105   ! max. call of simulator
      IF (mode.eq.6) info =  1106   ! max. precision
      IF (mode.eq.7) info = -1107   ! error matrix factorization
c
      RETURN
      END
