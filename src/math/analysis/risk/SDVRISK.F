c=======================================================================
c
c     subroutine SDVRISK                                     
c
c     Generates a portfolio history or future
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            njomax : number maximum of days
c            nasmax : maximum size of portfolio
c            njour  : estimation period
c            nasset : size of portfolio
c            tdata  : data of assets values during the estimation period
c                     matrix (njour*nasset)
c            porval : portfolio value at current date
c            wport  : weights in value of the portfolio (nasset)
c            test   : = 0, SDVRISK computes at the first date of tdata 
c                   : = 1, SDVRISK computes at the last  date of tdata 
c
c     OUTPUT :
c            histor : history of value of the portfolio (njour)
c            info   : =0 successful exit
c
c-----------------------------------------------------------------------
c
      subroutine SDVRISK ( njomax, nasmax, njour, nasset, tdata, porval,
     &                     wport, test, histor, info )
c
      implicit none
c
      integer njomax, nasmax, njour, nasset, test, info
      double precision porval
      double precision tdata(njomax,*)
      double precision wport(*), histor(*)
c
c     local variables
      integer i, j, k
      double precision somx, somz, aval
      double precision x(nasmax), z(nasmax)
c
c-----------------------------------------------------------------------
c
c     computing weights of portfolio in number of assets
      somx = 0.
c     if test=0 (first date), k=1; if test=1 (last date), k=njour 
      k = 1 + test*(njour-1)
      do i= 1,nasset
         aval = wport(i) / tdata(k,i)
         x(i) = aval
         somx =somx + aval
      end do
c
      do i= 1,nasset
         z(i) = x(i)/somx
      end do
c
c     computing history
      do j= 1,njour
         somz = 0.
c        if test=0 (first date), k=1 to njour
c        if test=1 (last date), k=njour to 1 
         k = j + test*(njour+1-2*j)
c
         do i= 1,nasset
            somz = somz + z(i) * tdata(k,i)
         end do
c
         histor(k) = somz * porval * somx
      end do
c
      info = 0
c
      return
      end
c
c=======================================================================
