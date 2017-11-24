c=======================================================================
c
c     subroutine pricing/PRvanilla                                      
c
c     This function computes analytic valuation and greeks for European 
c     style puts and calls (Vanilla) allows different daycount and interest 
c     rate conventions.
c
c     List of Results (Greeks)
c     ------------------------
c     0 :value
c     1 :spot delta
c     2 :gamma
c     3 :annual theta
c     4 :vega
c     5 :rho domestic
c     6 :rho foreign
c     7 :dstrike                 !derivative wrt strike
c     8 :d2strike                !second derivative wrt strike
c     9 :dduration               !derivative wrt to T(time to expiration) = - annual theta
c     10:leverage                !spot * spotdelta / value
c     11:vomma                   !second derivative wrt vol
c     12:retrieved_vola          !volatility implied by the option value, enter value for vol
c     13:forwarddelta            !derivative wrt forward
c     14:probdelta               !common probability factor in spot delta and forward delta, N(d1)
c     15:daily theta
c     16:forward
c     17:swaprate
c     18:d1
c     19:d2
c     20:disc                    !domestic discount factor
c     21:strike_of_delta     
c
c-----------------------------------------------------------------------
      SUBROUTINE PRvanilla(spot, strike, vol, rd, rf, tau, phi, 
     &                     convention, quotation, greek, van, info) 
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
c   convention : compounded convention                           integer
c                (0) linear for less than one year,
c                (1) annually compounded, 
c                (2) continuously compounded
c   quotation  : system of quotation                             integer
c                (0) domestic, 
c                (1) %domestic, 
c                (2) foreign, 
c                (3) %foreign
c   greek      : type of result                                  integer
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
c OUTPUT
c   van        : result of analytics                             double   
c   info       : diagnostic argument                            integer
c
c----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER phi, convention, quotation, greek, info
      DOUBLE PRECISION spot, strike, vol, rd, rf, tau, van
c
c     external functions
      DOUBLE PRECISION PRVAN
      EXTERNAL PRVAN
c
c----------------------------------------------------------------------
c
c     initialization 
      info = 0
      van  = 0
c
c     PRVAN (Analytics Vanilla), cf. pricing/PRutils.F    
      van = PRVAN(spot, strike, vol, rd ,rf ,tau ,phi ,convention,
     &             quotation, greek) 
c
      RETURN
      END   
 
 
