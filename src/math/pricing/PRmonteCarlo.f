c=======================================================================
c
c     subroutine PRmonteCarlo                                      
c
c     Monte Carlo Simulator for the geometric Brownian motion model
c     with 1) antithetic variables and 2) control variate
c 
c     Prices Vanilla and fixed strike Asian options
c
c-----------------------------------------------------------------------
      RECURSIVE DOUBLE PRECISION FUNCTION wy_mc(spot, strike, vol, rd, 
     & rf, tau, phi, convention, quotation, payofftype,
     & Nsim, Nsteps, greek) RESULT(f)      
c-----------------------------------------------------------------------
c
c INPUT
c   spot   : spot price (>0)                                      double  
c   strike : stike price (>0)                                     double
c   vol    : annual volatility (>=0)                              double
c   rd     : domestic risk free rate (>=0)                        double 
c            discounting is done as (1+r)^{-T}\n
c   rf     : foreign risk free rate or dividend rate (>=0)        double 
c            discounting is done as (1+rf)^{-T}\n 
c   tau    : time to maturity (in year)                           double
c   phi    : (1) Call, (-1) Put                                  integer
c   ea     : tree order (=1 or =2)                               integer 
c   convention : compounded convention                           integer
c                (0) linear for less than one year,
c                (1) annually compounded, 
c                (2) continuously compounded
c   quotation  : system of quotation                             integer
c                (0) domestic, 
c                (1) %domestic, 
c                (2) foreign, 
c                (3) %foreign
c   payofftype : (0) European put/call, 
c                (1) Goemetric fixed strike Asian put/call
c                (2) Arithmetic fixed strike Asian put/call
c   para1  : ?
c   para2  : ?
c   Nsim   : number of Monte Carlo simulations
c   Nsteps : number of steps ()
c   greek  : type of results cf. Results List 
c
c OUTPUT
c   mc         : result of Monte Carlo simulation                double   
c   info       : diagnostic argument                            integer
c
c----------------------------------------------------------------------
c
      INTEGER, PARAMETER :: Nstepsmax = 1000
      DOUBLE PRECISION, PARAMETER :: days = 365.0D0
      DOUBLE PRECISION :: spot, strike, vol, rd, rf, tau
      INTEGER :: phi, convention, quotation, payofftype, Nsim, Nsteps,
     &           greek, i, j
      DOUBLE PRECISION :: fc1, fc2, beta
      DOUBLE PRECISION :: my, volrt, delta_t, zeta,
     &                    f1, f2, geoavg1, geoavg2, ariavg1, ariavg2
      DOUBLE PRECISION, DIMENSION(0: Nstepsmax):: s1, s2
      
c
c     initial parameters      
      beta = 1.0D0          
      s1(0) = Log(spot) 
      s2(0) = Log(spot)
      delta_t = tau / days / Nsteps
      my = (rd - rf - 0.5D0 * vol * vol) * delta_t
      volrt = vol * Sqrt(delta_t)
      f1 = 0.D0
      f2 = 0.D0
      fc1 = 0.D0
      fc2 = 0.D0
      ! initialize random number generator
      CALL random_seed()

      DO j = 0, Nsim
        DO i = 1, Nsteps
            CALL random_number(zeta)
            
            !convert uniform zeta into standard normal zeta
            zeta = wy_ncinv(zeta)

            ! simulate path and antithetic path according to Euler scheme
            s1(i) = s1(i - 1) + my + volrt * zeta
            s2(i) = s2(i - 1) + my - volrt * zeta
        ENDDO
        SELECT CASE (payofftype)
        CASE (0) !European put/call
            f1 = f1 + wy_payoff(Exp(s1(Nsteps)), strike, phi)
            f2 = f2 + wy_payoff(Exp(s2(Nsteps)), strike, phi)
        CASE (1) !Goemetric fixed strike Asian put/call
          geoavg1 = 0.0D0
          geoavg2 = 0.0D0
          DO i = 0, Nsteps
                geoavg1 = geoavg1 + s1(i)
                geoavg2 = geoavg2 + s2(i)
          ENDDO
          SELECT CASE (greek)
          CASE (0) !value
            f1 = f1 + wy_payoff(Exp(geoavg1/(Nsteps+1)), strike, phi)
            f2 = f2 + wy_payoff(Exp(geoavg2/(Nsteps+1)), strike, phi)
          CASE (1) !delta
            IF (phi*geoavg1/(Nsteps+1) .gt. phi*log(strike)) THEN
                f1 = f1 + Exp(geoavg1/(Nsteps+1))/spot
            ENDIF
            IF (phi*geoavg2/(Nsteps+1) .gt. phi*log(strike)) THEN
                f2 = f2 + Exp(geoavg2/(Nsteps+1))/spot
            ENDIF
          END SELECT
        CASE (2) !Arithmetic fixed strike Asian put/call
          geoavg1 = 0.0D0
          geoavg2 = 0.0D0
          ariavg1 = 0.0D0
          ariavg2 = 0.0D0
          DO i = 0, Nsteps
            geoavg1 = geoavg1 + s1(i)
            geoavg2 = geoavg2 + s2(i)
            ariavg1 = ariavg1 + exp(s1(i))
            ariavg2 = ariavg2 + exp(s2(i))
          ENDDO
          fc1 = fc1 + wy_payoff(Exp(geoavg1/(Nsteps+1)), strike, phi)
          fc2 = fc2 + wy_payoff(Exp(geoavg2/(Nsteps+1)), strike, phi)
          f1  = f1 + wy_payoff(ariavg1/(Nsteps+1), strike, phi)
          f2  = f2 + wy_payoff(ariavg2/(Nsteps+1), strike, phi)
        END SELECT
      ENDDO
    
      ! get value of option with antithetic technique    
      f = 0.5D0 * (f1 + f2) / Nsim  * Exp(-rd * tau / days)
      IF (payofftype .eq. 2) THEN
        f = f + beta * (wy_Asiangeo(spot, strike, vol, rd, rf,
     &      tau / days, 0.0D0, spot, phi, convention, quotation, 0)
     &    - 0.5D0 * (fc1 + fc2) / Nsim  * Exp(-rd * tau / days))
      ENDIF
      END FUNCTION
