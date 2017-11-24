c=======================================================================
c
c     subroutine nqhess                                      
c
c     Quasi-Newton optimizer with BFGS method
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
      subroutine nqhess(n,imp,lp,iz,rz)
      implicit double precision(a-h,o-z)
      dimension iz(*),rz(*)
 1000 format(//)
 1001 format(34h nqhess   hessienne au point final)
 1002 format(9h   nqhess,i4,5d12.4,/,(9h   nqhess,4x,5d12.4))
      ni=2*n
      nw=n*(n+1)/2
      nw1=nw+n
      if(n.eq.1) go to 50
      nr=iz(ni+1)
      if(nr.eq.0) go to 20
      do 10 i=1,n
      if(iz(n+i).ne.0) go to 10
      nc=i
      call fajc1(n,nc,nr,rz(1),rz(nw+1),iz(1))
      if(nr.eq.0) go to 20
   10 continue
   20 n1=n-1
      do 40 i=1,n1
      j1=iz(i)
      if(j1.eq.i) go to 40
      ni=i
      nj=j1
      call f1qhes(n,ni,nj,nw,rz)
      call f1qhes(n,nj,ni,nw1,rz)
      call f2qhes(n,nj,nw,rz)
      call f2qhes(n,ni,nw1,rz)
      if(i.eq.n1) go to 50
      i1=i+1
      do 30 k=i1,n
      if(iz(k).ne.i) go to 30
      iz(k)=j1
      go to 40
   30 continue
   40 continue
   50 if(imp.le.0) return
      write(lp,1000)
      write(lp,1001)
      do 60 i=1,n
      write(lp,1002) i,(rz(i+(j-1)*(2*n-j)/2),j=1,i)
   60 continue
      return
      end
