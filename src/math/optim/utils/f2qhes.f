c=======================================================================
c
c     subroutine f2qhes                                      
c
c     Quasi-Newton optimizer with BFGS method
c
c-----------------------------------------------------------------------
c
      subroutine f2qhes(n,ni,nw,rz)      
      implicit double precision(a-h,o-z)
      dimension rz(*)
      nii=ni
      do 20 k=1,n
      rz(nii)=rz(nw+k)
      if(k.ge.ni) go to 10
      nii=nii+(n-k)
      go to 20
   10 nii=nii+1
   20 continue
      return
      end

