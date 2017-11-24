c=======================================================================
c
c     subroutine RETMU                                       
c
c     compute view on the diffusion parameter mu of an asset that
c     follows a GBM process with init value valinit and 
c     volatility parameter sigma from the view on the expected 
c     annualized return to maturity of the asset
c     method used = dichotomy
c
c-----------------------------------------------------------------------
c 
c     INPUT :
c            ret    : return to maturity                          double
c            nsimul : number of simulation                       integer
c            valini : initial value                               double
c            sigma  : volatility                                  double
c            matur  : maturity                                    double
c
c     OUTPUT :
c            mu     : diffusion parameter                         double
c            info   : = 0 successful exit                        integer
c                        
c     CALL   :
c        DICHOT  : DICHOTomy of a monotonous function
c                  on a gived interval
c                  for a gived y, find x than f(x) = y
c        MURET   : Computing diffusion parameter mu
c                  from return to maturity of the asset
c                  signature to be used by DICHOT
c                  (dichotomy of a monotonous function)
c
c-----------------------------------------------------------------------
c
      subroutine RETMU ( ret, nsimul, valini, sigma, matur,
     &                   mu, info)
c
      implicit none
c
      external MURET
c
      integer nsimul, info
      double precision ret, valini, sigma, matur, mu
c
      integer liudat, ldudat, maxitr, totitr
      integer iusdat(1)
      double precision xinf, xsup, epsabs
      double precision dusdat(3)
c
      parameter ( maxitr = 1000 )
      parameter ( xinf = -1., xsup = 1., epsabs = 1.0d-04 )
c
c----------------------------------------------------------------------
c
c     initialization
      liudat    = 1
      iusdat(1) = nsimul
      ldudat    = 3
      dusdat(1) = valini
      dusdat(2) = sigma
      dusdat(3) = matur
c
      call DICHOT ( MURET, ret, xinf, xsup, epsabs, maxitr,
     &              liudat, iusdat, ldudat, dusdat,
     &              mu, totitr, info )
c
      return
      end
