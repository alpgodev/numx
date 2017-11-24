c=======================================================================
c
c     subroutine SDVHST                                      
c
c     Generates a portfolio history or future with a risk-free asset.
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
c            wport  : weights in value of the portfolio (nasset+1)
c                     with risk-free asset as first component
c            seuil  : if weight of risk-free asset > seuil
c                     all the portfolio if risk-free
c            njan   : effective number of days in a year, interger
c            tasr   : yearly interest rate of risk-free asset, scalar
c            test   : = 0, SDVHST computes at the first date of tdata 
c                   : = 1, SDVHST computes at the last  date of tdata 
c
c     OUTPUT :
c            histor : history of value of the portfolio (njour)
c            info   : =0 successful exit
c
c     CALL   :
c            DAYRAT : computes a kind of daily 'rate' from the yearly
c                     interest rate of risk-free asset
c
c-----------------------------------------------------------------------
c
      subroutine SDVHST ( njomax, nasmax, njour, nasset, tdata, porval,
     &                    wport, seuil, njan, tasr, test, histor, info )
c
      implicit none
c
      integer njomax, nasmax, njour, nasset, njan, test, info
      double precision porval, seuil, tasr
      double precision tdata(njomax,*)
      double precision wport(*), histor(*)
c
c     local variables
      integer i, j, jdate
      double precision somx, somz, aval, taujou, mvpfin
      double precision x(nasmax), z(nasmax)
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     computing a kind of daily rate
      call DAYRAT ( tasr, njan, taujou, info )
      mvpfin = (njour-1)*taujou
c
c-----------------------------------------------------------------------
c
c     add risk asset
      if ( wport(1).lt.seuil) then
c
c     computing weights of portfolio in number of assets
      somx = 0.
      if ( test .eq. 0 ) then
         jdate = 1
      else
         jdate = njour
      end if
c
      do i= 1,nasset
         aval = wport(i+1) / tdata(jdate,i)
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
c
         somz = 0.
         do i= 1,nasset
            somz = somz + z(i) * tdata(j,i)
         end do
c
         histor(j) = porval*( somz*somx + wport(1) * (1.+(j-1)*taujou) 
     &                                            / ( 1 + test*mvpfin )
     &                      ) 
      end do
c
c     all with risk-free asset ( wport(1) = 1 )
      else
         do j= 1,njour
            histor(j) = porval * (1.+(j-1)*taujou )
     &                        / ( 1 + test*mvpfin )
         end do
      end if
      
      return
      end
c
c=======================================================================
