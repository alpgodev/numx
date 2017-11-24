c=======================================================================
c
c     subroutine  BFGSMVOLX (Expert Version)
c
c     Quasi-Newton optimizer with BFGS method
c
c-----------------------------------------------------------------------
      SUBROUTINE bfgsmvolx ( simul, simext, nvar, xvar,
     &                       dxmin, df1, epsabs, mode, totitr, totsim, 
     &                       binf, bsup, iwork, dwork,
     &                       funct, grad, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       simul  : entry point of the subroutine who computes
c                the function and the gradient
c       simext : entry point of an external subroutine
c       nvar   : number of variables                            integer
c       xvar   : initial value of the variables (nvar)           double
c       dxmin  : precision on variables stop criter (nvar)       double
c       df1    : criter decrease at the first iteration          double
c       epsabs : precision stop test                             double
c       mode   : initialization mode                            integer 
c       totitr : total of iterations                            integer
c       totsim : total of calls of SIMUL                        integer
c       binf   : inferior bounds on variables   vector(nvar)     double
c       bsup   : superior bounds on variables   vector(nvar)     double
c
c
c     WORKSPACE 
c       iwork  : 2*nvar+1                                        integer
c                + space for SIMUL
c       dwork  : nvar*(nvar+9)/2 
c                + space for SIMUL                                double
c                    
c     OUTPUT 
c       xvar   : variables value at the minimum,vector(nvar)      double
c       funct  : minimum of the function                          double
c       grad   : gradients value at the minimum,vector(nvar)      double
c       info   : = 0 successful exit                             integer
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
c     arguments i/o
      INTEGER nvar, mode, totitr, totsim, info
      DOUBLE PRECISION epsabs, df1, funct
      DOUBLE PRECISION xvar(*), binf(*), bsup(*), grad(*), dxmin(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)      
c
c     local variables
      INTEGER imp, ion2qn, indic, i, piwb, pdwb, piws, pdws
      REAL*4 rzs(1)
c
      integer iosort
      parameter ( iosort=6 )
c
      INTEGER minmax
      DOUBLE PRECISION normin, epsdxm, normg  
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
c      mode   = 1
c      totitr =  MAX(minmax,2*nvar)
c      totsim = totitr*3
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
      pdwb = 1
c     pdwb : pointer for QNBFGS who needs (nvar*(nvar+9)/2)
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
c      epsabs = epsbfg * nvar
c
c     computing df1 : estimation of the criter decrease
c                     at the first iteration
c                     must be nearly abs( f(xvar)-f(xopt) ) )
c      df1 = ABS(funct) / 10.
c
c     computing the precision on variables stop criter
c      normg = normin
c      DO i = 1,nvar
c         normg = MAX(ABS(grad(i)),normg)
c      ENDDO
c      DO i = 1,nvar
c        dwork(pdxmin+i-1) = ( df1* epsdxm ) / normg
c        dxmin(i) = ( df1* epsdxm ) / normg
c      ENDDO
c
c      open(unit=3,file='BFGSMVOLX_IN.txt',status='unknown')
c      write(3,*) "--- BEGIN BFGSMVOL ---"
c      write(3,*) "F=", funct
c      write(3,*) "gradF="
c      write(3,*) (grad(i), i=1,nvar)
c      write(3,*) "nb. iterations=", totitr
c      write(3,*) "nb. simulator=", totsim
c      write(3,*) "mode=", mode
c      write(3,*) "df1=", df1
c      write(3,*) "epsabs=", epsabs
c      write(3,*) "dxmin="
c      write(3,*) (dxmin(i), i=1,nvar)
c      write(3,*) "--- END BFGSMVOL ---"
c      close(unit=3)
c
c     optimization Quasi-Newton
      CALL qnbfgs ( simul, simext, nvar, xvar,
     &              funct, grad, dxmin, df1, epsabs,
     &              imp, ion2qn ,mode, totitr, totsim, binf, bsup,
     &              iwork(piwb), dwork(pdwb),
     &              iwork(piws), rzs, dwork(pdws) )
 
c      open(unit=4,file='BFGSMVOLX_OUT.txt',status='unknown')
c      write(4,*) "--- BEGIN BFGSMVOL ---"
c      write(4,*) "F=", funct
c      write(4,*) "gradF="
c      write(4,*) (grad(i), i=1,nvar)
c      write(4,*) "nb. iterations=", totitr
c      write(4,*) "nb. simulator=", totsim
c      write(4,*) "mode=", mode
c      write(4,*) "df1=", df1
c      write(4,*) "epsabs=", epsabs
c      write(4,*) "dxmin="
c      write(4,*) (dxmin(i), i=1,nvar)
c      write(4,*) "--- END BFGSMVOL ---"
c      close(unit=4)
c
c     error code
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
