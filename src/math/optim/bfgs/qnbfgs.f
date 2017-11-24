c=======================================================================
c
c     QNBFGS                                      
c
c     Quasi-Newton optimizer with BFGS method
c
c     Let f(.) a multidimentional function. 
c
c     Algorithm: x[i+1] = x[i] - inv(Hi)*Grad(f(x[i]))
c
c     Hi approximation of the Hesssian matrix of f(x)
c     limit (i->+oo) Hi = True Hessian matrix  
c
c-------------------------------------------------------------------------
      SUBROUTINE qnbfgs (simul, simext, n, x, f, g, dxmin, df1, epsabs,
     & imp, io, mode, iter, nsim, binf, bsup, iz, rz, izs, rzs, dzs)
c-------------------------------------------------------------------------
c
c     CALL:
c     n2qn1a, fmc11z, fajc1, fretc1, fmani1, fcomp1, fmlag1
c     nlis0, fcube, fuclid, fmc11a, fmc11b, fmc11e, nqhess, f1qhes, f2qhes
c
c-------------------------------------------------------------------------
c     
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      DIMENSION x(n), g(n), dxmin(n), bsup(n), binf(n)
      DIMENSION iz(*), rz(*), izs(*), dzs(*)
      REAL rzs(*)
c      
      EXTERNAL simul, simext
      DOUBLE PRECISION gpmopt
c      
      IF (imp .ne. 0) THEN
          nw = n*(9 + n)/2
          ni = 2*n + 1
          WRITE (io,1000) n,mode,iter,nsim,imp,df1,epsabs,ni,nw
      ENDIF
 1000 FORMAT (/" N2QN1 : point d'entree"/
     &    5x,"dimension du probleme (n):",i10/
     &    5x,"mode d'entree (mode):     ",i10/
     &    5x,"max iterations (iter):    ",i10/
     &    5x,"max simulation (nsim):    ",i10/
     &    5x,"niveau d'impression (imp):",i10/
     &    5x,"decroissancve attendue de f (df1):",1pd9.2/
     &    5x,"precision absolue (epsabs):       ",1pd9.2/
     &    /" Espace de travail:"/
     &    5x,"entiers:", i6/
     &    5x,"reels:  ", i6)
      IF (n .le. 0
     &  .or. (mode.eq.1 .and. df1.le.0.d0)
     &  .or. epsabs.lt.0.d0
     &  .or. imp.lt.0 .or. imp.gt.5
     &  .or. mode.lt.1 .or. mode.gt.4
     &  .or. iter.le.0
     &  .or. nsim.le.0) THEN
          WRITE (io,1001)
          mode=2
          RETURN
      ENDIF
 1001 FORMAT (/" >>> m2qn1: appel incoherent (n, df1, epsabs, imp,",
     &         " mode, iter ou nsim)"/)
      do i=1,n
          if (dxmin(i).le.0.d0
     &      .or. x(i).lt.binf(i)
     &      .or. x(i).gt.bsup(i)) then
              write (io,1002) i
              mode=2
              return
          endif
      enddo
c      
 1002 FORMAT (/" >>> n2qn1: borne(s) sur variable ",i4,
     &         " incoherente(s)"/)
c     
      nd     = 1 + (n*(n+1))/2
      nww    = nd + n
      nww1   = nww + n
      nga    = nww1 + n
      nindi  = 1
      nibloc = nindi+n
      ni     = nibloc+n
      s      = 0.d0
c
      DO 110 i=1,n
  110 s      = s + dxmin(i)*dxmin(i)
      epsabs = epsabs*dsqrt(s/dble(float(n)))
c
c     N2QN1A function call      
      CALL n2qn1a (simul,simext,n,x,f,g,dxmin,epsabs,df1,mode,
     & iter,nsim,imp,io,rz,rz(nd),rz(nww),rz(nww1),
     & rz(nga),binf,bsup,iz(nindi),iz(nibloc),iz(ni),
     & izs,rzs,dzs)
c     
      IF (imp.ge.2) WRITE(io,1003)
 1003 FORMAT (/1x,79("-"))
      IF (imp.ge.1) WRITE(io,1004) mode,iter,nsim,epsabs,gpmopt
 1004 FORMAT (" N2QN1: sortie en mode ",i2/
     &     5x,"nombre d'iterations   = ",i4/
     &     5x,"nombre de simulations = ",i4/
     &     5x,"|gradient projete| moyen       = ",1pd10.3/
     &     5x,"|gradient_dxmin projete| moyen = ",1pd10.3)
      IF (imp.ge.1) WRITE(io,1005) (i,iz(nibloc+i-1),i=1,n)
 1005 FORMAT (5x,"bornes",
     &    " (0: inactive, -1: binf active, +1: bsup active)",/
     &    (10x,"| ",5(i5,": ",i2," |")))
      RETURN
      END
