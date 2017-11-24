c=======================================================================
c
c     GENPRM                                                 
c 
c     Generate random permutation of array
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      SUBROUTINE genprm(array, larray)
c-----------------------------------------------------------------------      
c
c     INPUT : 
c     ARRAY <--> On output ARRAY is a random permutation of its
c                 value on input
c                         DOUBLE PRECISION ARRAY( LARRAY )
c     LARRAY <--> Length of ARRAY    integer
c
c-----------------------------------------------------------------------
c
      implicit none
c      
c     scalar arguments
      INTEGER larray
c     
c     array arguments
      DOUBLE PRECISION array(larray)
c     
c     local scalars
      INTEGER i, iwhich
      DOUBLE PRECISION elt, llarray
c     
c     external functions
      DOUBLE PRECISION ignuin
      EXTERNAL ignuin
c     
c     executable statements 
      llarray = dble(larray)
      DO i = 1,larray
         iwhich = int(ignuin(dble(i),llarray))
         elt = array(iwhich)
         array(iwhich) = array(i)
         array(i) = elt
      END DO
      END
