c=======================================================================
c
c     simulateur et fonction de trace de courbe pour RNLS
c
c     polynomial function f = an*x**n + .......... + a1*x + a0
c
c=======================================================================
c
c     SIMEXT
c
c     polynome: a(n)*x**n + .......... + a(1)*x + a(0)
c
c     Subroutine provided by the user to compute the value of
c     the functions vector(nfct) and the jacobian matrix(nfct*mx)
c     in function of xsol vector(mx)
c
c=======================================================================
      SUBROUTINE SIMEXT ( nfct, mxopt, xopt, spadim, grid,
     &                    liudat, iusdat, ldudat, dusdat,
     &                    funct, jacob, info )
c=======================================================================     
c
c     INPUT :
c            nfct   : number of functions                        integer
c            mxopt  : size of the vector solution                integer
c            xopt   : solution  vector(mxopt)                     double
c                     in the order a(0),a(1),........,a(n-1),a(n)
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            liudat : size of integers data provided by the user integer
c            iusdat : integers user data provided by the user
c                     vector(liudat)                             integer
c            ldudat : size of doubles data provided by the user  integer
c            dusdat : doubles user data provided by the user
c                     vector(ldudat)                              double
c
c     OUTPUT :
c            funct  : functions   vector(nfct)                    double
c            jacob  : jacobian    matrix(nfct*mx)                 double
c            info   : = 0 successful exit                        integer
c                     = else : it's your problem
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer nfct,mxopt,spadim,liudat,ldudat,info
      integer iusdat(*)
      double precision xopt(*),grid(*),dusdat(*),funct(*),jacob(nfct,*)
c
      integer i,j
      double precision z,r
c
c-----------------------------------------------------------------------
c
      info = 0
c
       DO i=1,nfct
         z = grid(i)
         CALL POLYN ( mxopt-1, xopt,  z, funct(i) )
         r = 1.
         DO j=1,mxopt
            jacob(i,j) = r
            r = r * z
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     SIMTRA
c
c     Computing RNLS files for graphic output for polynomial function
c
c=======================================================================
c
c     INPUT :
c            npc    : number of functions                        integer
c            mparam : number of parameters                       integer
c            pref   : reference parameters  vector(mparam)        double
c            prnls  : RNLS parameters  vector(mparam)             double
c            pnls   : NLS parameters  matrix(mparam,nmes)         double
c            npobs  : number of observations                     integer
c            pobs   : points of observation  vector(npobs)        double
c            obs    : mean observations  vector(npobs)            double
c            nmes   : number of measures                         integer
c            mesure : measures matrix(nmes,npobs)                 double
c
c-----------------------------------------------------------------------
c
      subroutine SIMTRA ( npc, mparam, pref, prnls, pnls, npobs, pobs,
     &                    obs, nmes, mesure )
c
      implicit none
c
      integer npc, mparam, npobs, nmes
      double precision pnls(mparam,*), mesure(nmes,*)
      double precision pref(*), prnls(*), pobs(*), obs(*)
c
c.......................................................................
c
      integer i, j, imes
      double precision pasc, z, debord, debut, fin
      double precision zc(npc), fc(npc,nmes), fref(npc), frnls(npc)
c
c-----------------------------------------------------------------------
c
      open(20,status='unknown',file='dim')
      open(21,status='unknown',file='obs')
      open(22,status='unknown',file='mod')
c
      write(20,*) npobs, nmes, npc
c
      debord = ( pobs(npobs) - pobs(1) ) / 100.
      debut =  pobs(1) - debord
      fin =  pobs(npobs) + debord
      pasc = ( fin - debut ) / ( npc - 1 )
c
      z = debut
      do i=1,npc
         zc(i) = z
         z = z + pasc
      end do
c
c     ecriture des points de mesure
c
      do i=1,npobs
         write(21,*)pobs(i),obs(i),(mesure(j,i),j=1,nmes)
      end do
c
c     calcul des nmes courbes
c
      do i=1,npc
         do imes=1,nmes
            call POLYN ( mparam-1, pnls(1,imes), zc(i), fc(i,imes) )
         end do
         call POLYN ( mparam-1, pref, zc(i), fref(i) )
         call POLYN ( mparam-1, prnls, zc(i), frnls(i) )
      end do
c
c     ecriture des courbes
c
      do i=1,npc
         write(22,*)zc(i),(fc(i,imes),imes=1,nmes),frnls(i),fref(i)
      end do
c
      close(20)
      close(21)
      close(22)
c
      return
      end
