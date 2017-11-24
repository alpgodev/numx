c=======================================================================
c
c     subroutine CORVOLM                                     
c
c     Empirical Vorrelation matrix, Volatility and Mean vector
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            ndate  : estimation period                          integer
c            nasset : size of portfolio                          integer
c            rasset : returns assets during the estimation period
c                     vectorized matrix (ndate*nasset)            double
c
c     WORKSPACE :
c            dwork  : vector( nasset*(ndate+nasset) )             double
c
c     OUTPUT :
c            rmean  : mean returns  vector(nasset)                double
c            volat  : volatilities  vector(nasset)                double
c            corr   : correlation matrix  vector(nasset*nasset)   double
c            info   : = 0 successful exit                        integer
c             
c     CALL   :
c        MCM     : computing the mean of each column of
c                  a vectorized full matrix
c        EMPCOV  : computes the empiric covariance matrix
c
c-----------------------------------------------------------------------
c
      subroutine corvolm ( ndate, nasset, rasset, dwork,
     &                     rmean, volat, corr, info)
c
      implicit none
c
      integer ndate, nasset, info
      double precision rasset(*), rmean(*), volat(*), corr(nasset,*)
      double precision dwork(*)
c
c     local variables
      integer prel, pcov, i, j, kij, kii
      double precision EPS, val
      parameter ( EPS = 1.E-30 )
c
c-----------------------------------------------------------------------
c
c     initialisations 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      prel = 1
c     prel   : pointer for relatives rents, so ( ndate*nasset ) more
      pcov = prel + ndate*nasset
c     pcov   : pointer for vectorized matrix
c              so (nasset*nasset) more
c
c     Total size of dwork array = nasset*(ndate+nasset) 
c
c-----------------------------------------------------------------------
c
c     computing means of rents
      call MCM ( ndate, nasset, rasset, rmean )
c
c     computing empiric covariance matrix
      call empcov ( ndate, nasset, rasset, rmean,
     &              dwork(prel), dwork(pcov), info )
      IF (info .LT. 0) RETURN
c
c     computing of square root of diagonal
      do i=1,nasset
         kii = pcov + (i-1)*(nasset+1)
         val = dwork(kii)
         if ( val .GT. EPS ) then
            volat(i)= SQRT(val)
         else
            info = -5
         end if
      end do
c
c     builds the empirical correlation matrix
      corr(1,1) = 1.
      do j=2,nasset
         do i=1,j
            kij = pcov + (i-1)*nasset + (j-1)
            corr(i,j) = dwork(kij)/( volat(i)*volat(j) )
            if ( i.ne.j )corr(j,i) = corr(i,j)
         end do
         corr(j,j) = 1.
      end do
c
      return
      end

