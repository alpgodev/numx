c=======================================================================
c
c     SIMSKEW                              
c
c     This function implements the simulator for the allocation with
c      skewness.
c
c-----------------------------------------------------------------------
      subroutine simskew ( indic, simext, ndual, ydual, funct, grad,
     &                     iwork, dwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       indic  : = 4, to compute function and gradients          integer
c       simext : entry point of an external subroutine
c                provided by the user (not used here)
c       ndual   : number of variables                            integer
c       ydual  : dual solution lambda (ndual)                     double
c
c     OUTPUT 
c       funct  : dual function value                              double
c       grad   : gradient value                                   double
c
c     WORKSPACE 
c
c       iwork : 5 + 3*n + neq + 2*nin                            integer
c
c       dwork : 3*(n*n) + 17*n + n*ntot + 2*ntot + 2*nin + 3     double
c           
c     CALL  
c       buildthgen, qp, skthgenopt
c
c-----------------------------------------------------------------------
c
      implicit none 
c
      external simext
c
c     i/o arguments
      integer indic, ndual
      double precision funct, ydual(*), grad(*)
c
c     workspaces      
      integer iwork(*)
      double precision dwork(*)
      integer pdQ, pdrho, pdw0, pdcinf, pdcsup, pdc, pdb
      integer pdcov, pdlin, pdwopt, pdlagr, pdw, pdqlin
      integer pdtemp, pdtemp2, pdtemp3 
      integer pdmu, pddelta, pddist
      
      integer pin, pinineq, pineq, piinfo, piwqua
c
c     local variables 
      integer info, i, n, neq, nin, ntot, nsq
      double precision temp, mu, delta, dist
c
c-----------------------------------------------------------------------
c
c     initializations
      
c
c     pointers for integer work space  : iwork
c     ---------------------------------------- 
      pin = 1
c     needs 1
c      so 1
      pinineq = pin + 1
c     needs 1
c      so 2
      pineq = pinineq + 1
c     needs 1 
c      so 3
      piinfo = pineq + 1
c     needs 1
c      so 4
      piwqua =  piinfo + 1
c     piwqua  : pointer for QUAPRO (3*n + neq + 2*nin + 1)
c             
c     Total so size of iwork array = 3*n + neq + 2*nin + 5

      info = 0
      n    = iwork(pin)
      nin  = iwork(pinineq)
      neq  = iwork(pineq)
      nsq = n*n
      ntot = nin + neq
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdmu = 1
c     needs 1
c     so 1
      pddelta = pdmu + 1  
c     needs 1
c     so 2
      pddist = pddelta + 1  
c     needs 1
c     so 3
      pdQ = pddist + 1 
c     needs n*n 
c     so n*n + 3
      pdqlin = pdQ + nsq
c     needs n 
c     so n*n + n + 3     
      pdrho = pdqlin + n
c     needs n
c     so n*n + 2*n + 3 
      pdw0 = pdrho + n
c     needs n
c     so n*n + 3*n + 3
      pdcinf = pdw0 + n
c     needs n
c     so n*n + 4*n + 3 
      pdcsup = pdcinf + n
c     needs n
c     so n*n + 5*n + 3 
      pdc = pdcsup + n 
c     needs n*ntot
c     so n*n + 5*n + n*ntot + 3
      pdb = pdc + n*ntot 
c     needs ntot
c     so n*n + 5*n + n*ntot + ntot + 3 
      pdtemp = pdb + ntot
c     needs n
c     so n*n + 6*n + n*ntot + ntot + 3 
      pdtemp2 = pdtemp + n
c     needs n
c     so n*n + 7*n + n*ntot + ntot + 3 
      pdtemp3 = pdtemp2 + n
c     needs n
c     so n*n + 8*n + n*ntot + ntot + 3 
      pdcov = pdtemp3 + n
c     needs n*n
c     so 2*(n*n) + 8*n + n*ntot + ntot + 3 
      pdlin = pdcov + n*n
c     needs n
c     so 2*(n*n) + 9*n + n*ntot + ntot + 3 
      pdwopt = pdlin + n
c     pdwopt : pointer for wopt n     
c     so 2*(n*n) + 10*n + n*ntot + ntot + 3  
      pdlagr = pdwopt + n
c     pdlagr : pointer for Lagrange mult. n+neq+nin
c     so 2*(n*n) + 11*n + n*ntot + 2*ntot + 3 
      pdw    = pdlagr + n + ntot
c     pdw    : pointer for QUAPRO (n*n + 6*n + 2*nin) 
c               pointer for SKTHGENOPT n                 
c     so 3*(n*n) + 17*n + n*ntot + 2*ntot + 2*nin + 3 
c
c-----------------------------------------------------------------------
c

c
c     Building the variables cov and lin for the optimization 
c     with quapro
      call buildthgen ( n, dwork(pdqlin), dwork(pdrho), dwork(pdw0),
     &                  dwork(pdQ), ndual, ydual, dwork(pdw),
     &                  dwork(pdcov), dwork(pdlin))
c      call imprv(6, ' cinf ', n, dwork(pdcinf))
c      call imprv(6, ' csup ', n, dwork(pdcsup))
      
      CALL qp( n, dwork(pdcov), dwork(pdlin), neq, nin, dwork(pdc),
     &         dwork(pdb), dwork(pdcinf), dwork(pdcsup),
     &         iwork(piwqua), dwork(pdw),
     &         dwork(pdlagr), dwork(pdwopt), iwork(piinfo))
      
c      call imprv(6, ' lagr ', n+neq, dwork(pdlagr))
c      call imprv(6, ' wopt ', n, dwork(pdwopt))
       call skthgenopt ( n, dwork(pdcov), dwork(pdlin), dwork(pdrho),
     &                   dwork(pdwopt), dwork(pdw0),
     &                   ndual, ydual, dwork(pddelta),
     &                   dwork(pdmu), dwork(pddist), dwork(pdw),
     &                   funct, grad)
     
c      write(6, *) 'funct ', funct
c      call imprv(6, 'grad ', ndual,  grad)
      
      return
      end
