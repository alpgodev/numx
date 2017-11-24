c=======================================================================
c
c     subroutine f1qhes                                      
c
c     Quasi-Newton optimizer with BFGS method
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
      subroutine f1qhes(n,ni,nj,nw,rz)
      implicit double precision(a-h,o-z)
      dimension rz(*)
      nii=ni
      nwi=nw+ni
      nwj=nw+nj
      do 20 k=1,n
      rz(nw+k)=rz(nii)
      if(k.ge.ni) go to 10
      nii=nii+(n-k)
      go to 20
   10 nii=nii+1
   20 continue
      rznw=rz(nwi)
      rz(nwi)=rz(nwj)
      rz(nwj)=rznw
      return
      end
