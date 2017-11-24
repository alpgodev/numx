c=======================================================================
c
c     subroutine PRutils                                      
c
c     Utils functions for pricing
c
c-----------------------------------------------------------------------
      recursive double precision function PRCRR(spot, strike, vol, 
     &    rd, rf, tau, phi, ea, convention, quotation, greek, N) 
     &    RESULT (f)
      PARAMETER (percent = 0.01, nmax = 1002)
      DOUBLE PRECISION :: spot, strike, vol, rd, rf, tau
      INTEGER ::ea, convention, quotation, greek, N, phi
      DOUBLE PRECISION, dimension(-nmax: nmax):: s
      DOUBLE PRECISION, dimension(-1: nmax):: V
      DOUBLE PRECISION :: dt, DF, u, p, q, relstrike, rforward, disc, 
     &                    valuefactor
      INTEGER :: i, j, k
     
      SELECT CASE (greek)      
      CASE (0) ! option value
        f = 0. ! initialize, then do two binomial trees, 
               ! one with N steps and one with N+1 steps
        DO k = 0, 1
            ! set parameters
            dt = tau / (N + k)
            u = exp(vol * Sqrt(dt))
 
            SELECT CASE (convention)
            CASE (0) ! linear for for mat less than one year, 
                     ! annually compounded for higher maturities
                IF (tau .le. 1.) THEN
                    rforward = (1. + rd * tau * 365. / 360.) / 
     &                         (1. + rf * tau * 365. / 360.)
                    disc = 1. / (1. + rd * tau * 365. / 360.)
                ELSE
                    rforward = exp((log(1. + rd) - log(1. + rf)) * tau)
                    disc     = exp(-log(1 + rd) * tau)
                ENDIF
            CASE (1) ! annually compounded
                rforward = exp((log(1. + rd) - log(1. + rf)) * tau)
                disc     = exp(-log(1. + rd) * tau)
            CASE (2) ! continuously compounded
                rforward = exp((rd - rf) * tau)
                disc     = exp(-rd * tau)
            ENDSELECT
 
            DF = exp(log(disc) / (N+k))
            p  = (exp(log(rforward) /(N+k)) - 1. / u) / (u - 1. / u)*DF
            q = DF - p
            relstrike = strike / spot
            ! build the stock price tree, assuming spot = 1:
            DO j = -(N + k), N + k
                s(j) = u**j
            ENDDO
            ! set final time value:
            DO j = 0, N + k
                V(j) = wy_payoff((s(2 * j - N - k), relstrike, phi)
            ENDDO
            ! compute present value by stepping backwards in the value tree:
            DO i = N + k, 1,-1
                DO j = 0, i - 1
                    V(j) = (p * V(j + 1) + q * V(j))
                    if (ea .eq. 2) then
                    V(j) = max(V(j), wy_payoff((s(2 * j - i + 1)),
     &                                          relstrike, phi))
                    end if
                ENDDO
            ENDDO
            f = f + V(0)
        ENDDO
        f = f * spot / 2.
        
      CASE (1) ! delta
      f = (PRCRR(spot*(1.+percent) , strike, vol, rd, rf, tau, phi,ea,
     &           convention ,0 , 0, N)
     &  - PRCRR(spot*(1.-percent) , strike, vol, rd, rf, tau, phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * spot)
      valuefactor = 1 / spot
      CASE (2) !gamma
      f = (PRCRR(spot*(1.+percent) , strike, vol, rd, rf, tau, phi, ea,
     &           convention ,0 , 0, N)
     &    + PRCRR(spot*(1.-percent) , strike, vol, rd, rf, tau, phi, ea,
     &           convention ,0 , 0, N)
     &   - 2.*PRCRR(spot, strike, vol, rd, rf, tau, phi, ea, convention,
     &              0 , 0, N)) / (percent*percent*spot*spot)
      CASE (3) !theta
      f = (PRCRR(spot, strike, vol, rd, rf, tau*(1.-percent) , phi, ea,
     &           convention ,0 , 0, N)
     &  - PRCRR(spot, strike, vol, rd, rf, tau*(1.+percent) , phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * tau)
      CASE (4) !vega
      f = (PRCRR(spot, strike, vol*(1.+percent) , rd, rf, tau, phi, ea,
     &           convention ,0 , 0, N) 
     &  - PRCRR(spot, strike, vol*(1.-percent) , rd, rf, tau, phi, ea,
     &          convention ,0 , 0, N)) / (200. * percent * vol)
      CASE (5) !rhod
      f = (PRCRR(spot, strike, vol, rd*(1.+percent) , rf, tau, phi, ea,
     &           convention ,0 , 0, N) 
     &  - PRCRR(spot, strike, vol, rd*(1.-percent) , rf, tau, phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * rd)
      CASE (6) !rhof
      f = (PRCRR(spot, strike, vol, rd, rf*(1.+percent) , tau, phi, ea,
     &           convention ,0 , 0, N) 
     &  - PRCRR(spot, strike, vol, rd, rf*(1.-percent) , tau, phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * rf)
      CASE (7) !dstrike
      f = (PRCRR(spot, strike*(1.+percent) , vol, rd, rf, tau, phi, ea,
     &           convention ,0 , 0, N) 
     &  - PRCRR(spot, strike*(1.-percent) , vol, rd, rf, tau, phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * strike)
      CASE (8) !d2strike
      f = (PRCRR(spot, strike*(1.+percent) , vol, rd, rf, tau, phi, ea,
     &           convention ,0 , 0, N) 
     &  + PRCRR(spot, strike*(1.-percent) , vol, rd, rf, tau, phi, ea, 
     &          convention ,0 , 0, N) 
     &  - 2.*PRCRR(spot, strike, vol, rd, rf, tau, phi, ea, convention,
     &             0 , 0, N)) / (percent*percent*strike*strike)
      CASE (9) !dduration
      f = (PRCRR(spot, strike, vol, rd, rf, tau*(1.+percent) , phi, ea,
     &           convention ,0 , 0, N) 
     &  - PRCRR(spot, strike, vol, rd, rf, tau*(1.-percent) , phi, ea,
     &          convention ,0 , 0, N)) / (2. * percent * tau)
      CASE (10)!leverage
      f = spot * (PRCRR(spot*(1.+percent) , strike, vol, rd, rf, tau,
     &                  phi, ea, convention ,0 , 0, N) 
     &  - PRCRR(spot*(1.-percent) , strike, vol, rd, rf, tau, phi, ea, 
     &          convention ,0 , 0, N)) / (2. * percent * spot) 
     &  / PRCRR(spot, strike, vol, rd, rf,tau,phi,ea,convention,0, 0, N)
      CASE (11)! vomma
      f = (PRCRR(spot, strike, vol*(1.+3.*percent) , rd, rf, tau, phi, 
     &           ea, convention ,0 , 0, N) 
     &  + PRCRR(spot, strike, vol*(1.-3.*percent) , rd, rf, tau, phi, 
     &          ea, convention ,0 , 0, N) 
     &  - 2.*PRCRR(spot, strike, vol, rd, rf, tau, phi, ea, convention,
     &             0 , 0, N)) / (9.*percent*percent*vol*vol)
      ENDSELECT
 
!quotation adjustment:
      SELECT CASE (greek)
      CASE (0, 3, 4, 5, 6, 9)
        SELECT CASE (quotation)
        CASE (0) !domestic
        CASE (1) !%domestic
            f = f / strike * 100.
        CASE (2) !foreign
            f = f / strike / spot
        CASE (3) !%foreign
            f = f / spot * 100.
        END SELECT
      CASE (1, 7)
        SELECT CASE (quotation)
        CASE (0) !plain derivative
        CASE (1) !plain derivative in %
            f = f * 100.
        CASE (2) !spot * plain derivative of %foreign quoted value
            f = f - PRCRR(spot, strike, vol, rd, rf, tau, phi, ea, 
     &                    convention, 0, 0, N) * valuefactor
        CASE (3) !spot * plain derivative of %foreign quoted value in%
            f = 100. * (f - PRCRR(spot, strike, vol, rd, rf, tau, phi,
     &                           ea, convention, 0, 0, N) * valuefactor)
        END SELECT
      CASE (2)
        SELECT CASE (quotation)
        CASE (0) !plain derivative
        CASE (1) !change of delta when spot changes by 1%
            f = f * spot
        END SELECT
      CASE (8)
        SELECT CASE (quotation)
        CASE (0) !plain derivative
        CASE (1) !change of delta when strike changes by 1%
            f = f * strike
        END SELECT
      END SELECT
 
      END FUNCTION
 
!---------------------------------------------------------------------------------------------------------------
!     Payoff function 
!
!     spot   : spot price
!     strike : strike price 
!     phi    : (1) call payoff, (-1) put payoff
!---------------------------------------------------------------------------------------------------------------      
      DOUBLE PRECISION FUNCTION wy_payoff(spot, strike, phi)

      DOUBLE PRECISION spot, strike
      INTEGER phi
      
      ! payoff evaluation: max(0, phi*(S-K))  
      wy_payoff = max(0.D00, phi * (spot - strike))
      
      END FUNCTION

!-----------------------------------------------------------------------
!     Analytics Vanilla
!-----------------------------------------------------------------------
      recursive double precision Function PRVAN(spot , strike ,vol,
     & rd , rf , tau , phi , convention, quotation , greek) 
     & result(f) 
c-----------------------------------------------------------------------
c     Inputs
c       spot       : spot price
c       strike     : stike price
c       vol        : annual volatility
c       rd         : domestic risk free rate: discounting is done as (1+r)^{-T}\n
c       rf         : foreign risk free rate or dividend rate: discounting is done as (1+rf)^{-T}\n 
c       tau        : time to maturity (in year)
c       phi        : (1) Call, (-1) Put
c       convention : (0), (1) annually compounded, (2) continuously compounded
c       quotation  : (0) domestic, (1) %domestic, (2) foreign, (3) %foreign
c                    ...
c       greek      : cf. results list
c     (0) :value
c     (1) :spot delta, first derivative of the option value with respect to the underlying price
c     (2) :gamma, second derivative of the option value with respect to the underlying price
c     (3) :annual theta, derivative of the option value with respect to the amount of time to expiry
c     (4) :vega, derivative of the option value with respect to the volatility of the underlying
c     (5) :rho domestic, derivative of the option value with respect to the risk free rate
c     (6) :rho foreign
c     (7) :dstrike, derivative wrt strike
c     (8) :d2strike, second derivative wrt strike
c     (9) :dduration, derivative wrt to T(time to expiration) = - annual theta
c     (10):leverage, spot * spotdelta / value
c     (11):vomma/volga, second derivative wrt vol (second derivative of the option value with respect to the volatility of the underlying)
c     (12):retrieved_vola, volatility implied by the option value, enter value for vol
c     (13):forwarddelta, derivative wrt forward
c     (14):probdelta, common probability factor in spot delta and forward delta, N(d1)
c     (15):daily theta
c     (16):forward
c     (17):swaprate
c     (18):d1
c     (19):d2
c     (20):disc, domestic discount factor
c     (21):strike_of_delta 
c
c----------------------------------------------------------------------
    
      double precision :: spot, strike, vol, rd, rf, tau
      integer :: phi, convention, quotation, greek
      double precision :: d1, d2, forward, disc, convrd, convrf, 
     &                    dconvrd, dconvrf, valuefactor
 
      Select Case (convention)
      Case (0) 
      !linear for for mat less than one year, annually compounded for higher maturities
      If (tau .lt. 1.) Then
        forward = spot * (1. + rd * tau * 365. / 360.) 
     &                 / (1. + rf * tau * 365. / 360.)
        disc = 1. / (1. + rd * tau * 365. / 360.)
        convrd  = (rd * 365. / 360.) / (1. + rd * tau * 365. / 360.)
        convrf  = (rf * 365. / 360.) / (1. + rf * tau * 365. / 360.)
        dconvrd = (365. / 360.) / (1. + rd * tau * 365. / 360.)
        dconvrf = (365. / 360.) / (1. + rf * tau * 365. / 360.)
      Else
        forward = spot * Exp((Log(1. + rd) - Log(1. + rf)) * tau)
        disc = Exp(-Log(1. + rd) * tau)
        convrd = Log(1. + rd)
        convrf = Log(1. + rf)
        dconvrd = 1. / (1. + rd)
        dconvrf = 1. / (1. + rf)
      End If
      Case (1) !annually compounded
        forward = spot * Exp((Log(1. + rd) - Log(1. + rf)) * tau)
        disc = Exp(-Log(1. + rd) * tau)
        convrd = Log(1. + rd)
        convrf = Log(1. + rf)
        dconvrd = 1. / (1. + rd)
        dconvrf = 1. / (1. + rf)
      Case (2) !continuously compounded
        forward = spot * Exp((rd - rf) * tau)
        disc = Exp(-rd * tau)
        convrd = rd
        convrf = rf
        dconvrd = 1.
        dconvrf = 1.
      End Select
      
      If (greek .ne. 21) Then
        d1 = (Log(forward / strike) + vol * vol * 0.5 * tau) 
     &       / (vol * Sqrt(tau))
        d2 = d1 - vol * Sqrt(tau)
      End If
    
      Select Case (greek)
      Case (0) !value
        f = phi * disc * (forward * wy_nc(phi * d1))
     &      - strike * wy_nc(phi * d2))
      Case (1) !spotdelta
        f = phi * disc * forward / spot * wy_nc(phi * d1)
        valuefactor = 1. / spot
      Case (2) !gamma
        f = disc*forward / spot * wy_ndf(d1) / (vol * Sqrt(tau) * spot)
        Select Case (quotation)
        Case (0) !plain derivative
        Case (1) 
        !change in delta when spot changes by 1% in % (trader's Gamma)
            f = f * spot
        End Select
      Case (3) !annual theta
        f = -disc * (wy_ndf(d1) * vol * forward / (2. * Sqrt(tau))
     &      + phi * (convrd * strike * wy_nc(phi * d2) 
     &      - convrf * forward * wy_nc(phi * d1)))
      Case (4) !vega
        f = disc * forward * wy_ndf(d1) * Sqrt(tau) / 100.
      Case (5) !domestic rho
        f = phi * strike * tau * disc * wy_nc(phi * d2) * dconvrd
      Case (6) !foreign rho
        f = -phi * forward * tau * disc * wy_nc(phi * d1) * dconvrf
      Case (7) !dstrike
        f = -phi * disc * wy_nc(phi * d2)
        valuefactor = 0.
      Case (8) !d2strike
        f = disc * wy_ndf(d2) / (strike * vol * Sqrt(tau))
      Case (9) !dduration
        f = disc * (wy_ndf(d1) * vol * forward / (2. * Sqrt(tau)) 
     &    + phi * (convrd * strike * wy_nc(phi * d2) 
     &    - convrf * forward * wy_nc(phi * d1)))
      Case (10) !leverage
        f = forward * wy_nc(phi * d1) / (forward * wy_nc(phi * d1) 
     &    - strike * wy_nc(phi * d2))
      Case (11) !vomma
        f = disc * forward * wy_ndf(d1) * Sqrt(tau) * d1 * d2 / vol
      Case (12) !retrieved_vola
        f = wy_find_vol(spot, strike, vol, rd, rf, tau, phi, convention,
     &                  quotation, 0.1D0)
        Select Case (quotation)
        Case (0) !plain volatility
        Case (1) !volatility in %
            f = f * 100.
        End Select
      Case (13) !forwarddelta
        f = phi * wy_nc(phi * d1) * disc
        valuefactor = 1. / forward
      Case (14) !probdelta
        f = phi * wy_nc(phi * d1)
        valuefactor = 1. / (forward * disc)
      Case (15) !daily theta
        f = (-disc * (wy_ndf(d1) * vol * forward / (2. * Sqrt(tau)) 
     &    + phi * (convrd * strike * wy_nc(phi * d2) 
     &    - convrf * forward * wy_nc(phi * d1)))) / 365.
      Case (16) !forward
        f = forward
      Case (17) !swaprate
        f = (forward - spot) * 10000.
      Case (18) !d1
        f = d1
      Case (19) !d2
        f = d2
      Case (20) !domestic discount factor
        f = disc
      Case (21) !strike_of_delta (spot delta is given, strike is wanted)
        Select Case (quotation)
        Case (0) !plain derivative
            f = forward * Exp(-phi * wy_ncinv(phi * strike / disc 
     &        / forward*spot)*vol * Sqrt(tau) + 0.5D0 * vol * vol * tau)
        Case (1) !plain derivative in %
            f = forward * Exp(-phi * wy_ncinv(phi * strike / 100.D0/ 
     &          disc/forward*spot)*vol*Sqrt(tau) + 0.5D0*vol*vol*tau)
        Case (2) !spot * plain derivative of %foreign quoted value
c            f = wy_strike_of_pfdelta(spot, strike, vol, disc, forward,
c     &                               tau, phi)
        Case (3) !spot * plain derivative of %foreign quoted value in%
c            f = wy_strike_of_pfdelta(spot, strike / 100.D0, vol, disc,
c     &                               forward, tau, phi)
        End Select
      End Select
    
      Select Case (greek)
      Case (0, 3, 4, 5, 6, 9, 15)
        Select Case (quotation)
        Case (0) !domestic
        Case (1) !%domestic
            f = f / strike * 100.
        Case (2) !foreign
            f = f / strike / spot
        Case (3) !%foreign
            f = f / spot * 100.
        End Select
      Case (1, 7, 13, 14)
        Select Case (quotation)
        Case (0) !plain derivative
        Case (1) !plain derivative in %
            f = f * 100.
        Case (2) !spot * plain derivative of %foreign quoted value
            f = f - phi * disc * (forward * wy_nc(phi * d1) - strike 
     &                            * wy_nc(phi * d2)) * valuefactor
        Case (3) !spot * plain derivative of %foreign quoted value in%
            f = 100. * (f - phi * disc * (forward * wy_nc(phi * d1) 
     &                    - strike * wy_nc(phi * d2)) * valuefactor)
        End Select
      End Select
      End Function

!------------------------------------------------------------------------------------------------------------------
! find_vol finds the volatility given a value of a vanilla option based on Newton's method
!------------------------------------------------------------------------------------------------------------------

      RECURSIVE DOUBLE PRECISION FUNCTION wy_find_vol(spot, strike, vol,
     & rd, rf, tau, phi, convention, quotation, startvol) RESULT(f)
!
!
      DOUBLE PRECISION :: spot, strike, vol, rd, rf, tau, startvol
      INTEGER :: phi, convention, quotation
      DOUBLE PRECISION :: func, dfunc
      INTEGER :: counter
 !
 !
      SELECT CASE (quotation) !invert the quotation:
      CASE (0) !domestic
      CASE (1) !%domestic
        vol = vol * strike / 100.D0
      CASE (2) !foreign
        vol = vol * strike * spot
      CASE (3) !%foreign
        vol = vol * spot / 100.D0
      END SELECT

      !start Newton's method:
      f = startvol
      func = PRVAN(spot, strike, startvol, rd, rf, tau, phi, 
     &             convention, 0, 0) - vol
      dfunc = PRVAN(spot, strike, startvol, rd, rf, tau, phi, 
     &              convention, 0, 4) * 100.D0
      counter = 0
      !Newton's loop: x(n+1) = x(n) - F[x(n)]/F'[x(n)]
      DO WHILE (Abs(func / dfunc) .ge. 1.0D-7 .and. counter .le. 20)
        f = f - func / dfunc
        func = PRVAN(spot, strike, f, rd, rf, tau, phi, convention,
     &               0, 0) - vol
        dfunc = PRVAN(spot, strike, f, rd, rf, tau, phi,convention,
     &                0, 4) * 100.D0
        counter = counter + 1
      ENDDO
      END FUNCTION 
      
!---------------------------------------------------------------------------------------------------------------
! nc returns the cumulative distribution function of a standard normal random variable
!---------------------------------------------------------------------------------------------------------------
      recursive double precision function wy_nc(x) result (f)

      double precision :: x
      double precision, dimension(1: 5):: a
      If (x .lt. -7.D0) Then
        f = wy_ndf(x) / sqrt(1.D0 + x * x)
      ElseIf (x .gt. 7.D0) Then
        f = 1.D0 - wy_nc(-x)
      Else
        f    = 0.2316419D0
        a(1) = 0.31938153D0
        a(2) = -0.356563782D0
        a(3) = 1.781477937D0
        a(4) = -1.821255978D0
        a(5) = 1.330274429D0
        f = 1.D0 / (1.D0 + f * Abs(x))
        f = 1.D0 - wy_ndf(x) * (a(1) * f + a(2) * f ** 2.D0 + a(3) 
     &      * f ** 3.D0 + a(4) * f ** 4.D0 + a(5) * f ** 5.D0)
        If (x .le. 0.D0) Then
            f = 1.D0 - f
        End if
      End If
      End Function

!---------------------------------------------------------------------------------------------------------------
! ndf returns the density function of a standard normal random variable
      double precision function wy_ndf(x)
      double precision :: x
      
      wy_ndf = 0.398942280401433D0 * exp(-x * x * 0.5D0) 
      !0.398942280401433 = 1/sqareroot(2*pi)
      end function
!---------------------------------------------------------------------------------------------------------------
 
 
!---------------------------------------------------------------------------------------------------------------
! ncinv returns the inverse of a standard normal distribution function (nc in this library)
!---------------------------------------------------------------------------------------------------------------
      double precision function wy_ncinv(x)

      double precision :: x
      double precision :: func, dfunc
      if (x .le. 0.D0) then
        wy_ncinv = -1.5D1
      elseif (x .ge. 1.D0) then
        wy_ncinv = 1.5D1
      else
        !start Newton's method:
        wy_ncinv = 0.D0
        func = 0.5D0 - x
        dfunc = 0.398942280401433D0
        Do while (Abs(func / dfunc) .ge. 1.D-7)
                wy_ncinv = wy_ncinv - func / dfunc
                func = wy_nc(wy_ncinv) - x
                dfunc = wy_ndf(wy_ncinv)
        end do
      endif
      end function
!---------------------------------------------------------------------------------------------------------------
 

!CARDANO produces the roots of the equation x³+ax²+bx+c=0.
!output: 0 : discriminant D which classifies the solution:
!        D>0: one real y(1) and two complex conjugate solutions
!             y(2)+y(3)i and y(2)-y(3)i
!        D=0: three real solutions including one double solution (y(3)=0)
!        D<0: three distinct real solutions y(1)<y(2)<y(3)
!output: 1,2,3 : y(output)
!reference: Harris/Stocker, Handbook of Mathematics and Computational Science,
!           Springer 1998, Chapter 2.5 Cubic Equations
!-----------------------------------------------------------------------------
      double precision function wy_cardano(a,b,c,output)
      
      DOUBLE PRECISION EPSILON
      parameter (epsilon = 0.000001)
      double precision a,b,c
      double precision p, q, D, phi, u, v
      double precision x(3)
      double precision y(3)
      double precision radicand
      integer output
      
      p=(3.*b-a*a)/3.
      q=c+2.*a*a*a/27.-a*b/3.
      D=(P/3.)**3.+q*q/4.
      
      !the irreducible case has a trigonometric form:
      if (D .LT. 0.) then
        phi=acos(-q/(2*sqrt((abs(p)/3)**3)))
        x(1) = -a/3.+2*sqrt(abs(p)/3)*cos(phi/3.)
        x(2) = -a/3.-2*sqrt(abs(p)/3)*cos((phi-3.141592654)/3.)
        x(3) = -a/3.-2*sqrt(abs(p)/3)*cos((phi+3.141592654)/3.)

        !sort solutions according to order:
        y(1) = min(x(1),x(2),x(3))
        y(3) = max(x(1),x(2),x(3))
        y(2) = x(1)+x(2)+x(3)-y(1)-y(3)
      else
        radicand =abs(-q/2.+sqrt(D))
        if (radicand .ge. epsilon) then
            u = sign(1.,-q/2.+sqrt(D))*exp(log(radicand)/3.)
        else
            u = 0
        end if
            radicand =abs(-q/2.-sqrt(D))
            if (radicand .ge. epsilon) then
                v=sign(1.,-q/2.-sqrt(D))*exp(log(radicand)/3.)
            else
                v = 0
            end if
            y(1) = -a/3.+u+v
            y(2) = -a/3.-(u+v)/2.
            y(3) = sqrt(3.)*(u-v)/2.
        end if
        select case (output)
        case (0)
                wy_cardano = D
        case (1:3)
                wy_cardano = y(output)
        end select
      end function

!------------------------------------------------------------------------------------------------------------
 

!------------------------------------------------------------------------------------------------------------
! normal random generator
!------------------------------------------------------------------------------------------------------------
      double precision function wy_normrv()

      double precision :: x

      call random_number(x)
      wy_normrv = wy_ncinv(x)
      
      end function

!------------------------------------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------------------------------------
! Asiangeo computes value and greek of fixed strike European puts and calls on the continuously sampled
! geometric average taken between -s and T (-backtau and tau), check Wilmott:Derivatives
!------------------------------------------------------------------------------------------------------------
 
      recursive double precision function wy_Asiangeo(spot, strike, vol,
     & rd, rf, tau, backtau, currentavg, phi, convention, quotation,
     & greek) result(f)

      double precision :: spot, strike, vol, rd, rf, tau, backtau, currentavg
      integer :: phi, convention, quotation, greek
      double precision :: d1, d2, forward, disc, convrd, convrf,
     &                   dconvrd, dconvrf, valuefactor, geofactor, alpha

      alpha = tau / (tau + backtau)
      Select Case (convention)
      Case (0) 
      !linear for for mat less than one year, annually compounded for higher maturities
      If (tau .le. 1.D0) Then
        forward = spot * (1.0D0 + rd * tau * 365.0D0 / 360.0D0) 
     &                 / (1.0D0 + rf * tau * 365.0D0 / 360.0D0)
        disc = 1.0D0 / (1.0D0 + rd * tau * 365 / 360)
        convrd = (rd * 365.0D0 / 360.0D0) 
     &         / (1.0D0 + rd * tau * 365.0D0 / 360.0D0)
        convrf = (rf * 365.0D0 / 360.0D0) 
     &         / (1.0D0 + rf * tau * 365.0D0 / 360.0D0)
        dconvrd = (365.0D0 / 360.0D0) 
     &          / (1.0D0 + rd * tau * 365.0D0 / 360.0D0)
        dconvrf = (365.0D0 / 360.0D0) 
     &          / (1.0D0 + rf * tau * 365.0D0 / 360.0D0)
        geofactor = Exp(-0.5D0 * tau * alpha * (Log(1.0D0 + rd) 
     &            - Log(1.0D0 + rf) + vol * vol 
     &            / 2.0D0 * (1.0D0 - 2.0D0 * alpha / 3.0D0)))
      Else
        forward    = spot * Exp((Log(1.0D0 + rd) - Log(1.0D0 + rf))*tau)
        disc       = Exp(-Log(1.0D0 + rd) * tau)
        convrd     = Log(1.0D0 + rd)
        convrf     = Log(1.0D0 + rf)
        dconvrd    = 1.0D0 / (1.0D0 + rd)
        dconvrf    = 1.0D0 / (1.0D0 + rf)
        geofactor  = Exp(-0.5D0 * tau * alpha * (Log(1.0D0 + rd) 
     &             - Log(1.0D0 + rf) + vol * vol / 2.0D0 
     &                * (1.0D0 - 2.0D0 * alpha / 3.0D0)))
      End If
      Case (1) !annually compounded
        forward = spot * Exp((Log(1.0D0 + rd) - Log(1.0D0 + rf)) * tau)
        disc    = Exp(-Log(1.0D0 + rd) * tau)
        convrd  = Log(1.0D0 + rd)
        convrf  = Log(1.0D0 + rf)
        dconvrd = 1.0D0 / (1.0D0 + rd)
        dconvrf = 1.0D0 / (1.0D0 + rf)
        geofactor = Exp(-0.5D0 * tau * alpha * (Log(1.0D0 + rd) 
     &            - Log(1.0D0 + rf) + vol * vol / 2.0D0 
     &               * (1.0D0 - 2.0D0 * alpha / 3.0D0)))
      Case (2) !continuously compounded
        forward   = spot * Exp((rd - rf) * tau)
        disc      = Exp(-rd * tau)
        convrd    = rd
        convrf    = rf
        dconvrd   = 1.0D0
        dconvrf   = 1.0D0
        geofactor = Exp(-0.5D0*tau*alpha*
     &                 (rd-rf+vol*vol/2.0D0*(1.0D0-2.0D0*alpha/3.0D0)))
      End Select

      Select Case (greek)
      Case (0) !value
        f = disc ** (1.0D0 - alpha) 
     &    * PRVAN(geofactor*spot**alpha*currentavg**(1.0D0-alpha),
     &      strike, alpha * vol / Sqrt(3.0D0), alpha * rd, alpha * rf,
     &      tau, phi, convention, 0, 0)
      Case (1) !spotdelta
        f= disc ** (1.0D0 - alpha) 
     &   * PRVAN(geofactor*spot**alpha*currentavg**(1.0D0-alpha), 
     &     strike, alpha * vol / Sqrt(3.0D0), alpha * rd, alpha * rf,
     &     tau, phi, convention, 0, 1)
     &     *geofactor*alpha*(currentavg/spot)**(1.0D0 - alpha)

      End Select
      End Function
