c=======================================================================
c
c     subroutine PRbinomial                                      
c
c     This function produces values and greeks of European and American 
c     puts and calls following the Cox, Ross, Rubinstein binomial model.
c
c-----------------------------------------------------------------------
      SUBROUTINE PRbinomial(spot, strike, vol, rd, rf, tau, phi, ea, 
     &                      convention, quotation, greek, N, crr, info) 
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
c   greek      : type of result                                  integer
c                (0) option pricing, 
c                (1) spot delta, 
c                (2) gamma, 
c                (3) annual theta, 
c                (4) vega, 
c                (5) rho domestic, 
c                (6) rho foreign, 
c                (7) derivative wrt strike, 
c                (8) second derivative wrt strike, ...
c                (9) duration or derivative wrt to T(time to expiration) = - annual theta, ...
c                (10) leverage = spot * spotdelta / value, 
c                (11) vomma or second derivative wrt vol
c   N          : step of the binomial trees                     integer    
c
c OUTPUT
c   crr        : result of CRR binomial model                    double   
c   info       : diagnostic argument                            integer
c
c----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER phi, ea, convention, quotation, greek, N, info
      DOUBLE PRECISION spot, strike, vol, rd, rf, tau, crr
c
c     external functions
      DOUBLE PRECISION PRCRR
      EXTERNAL PRCRR
c     
c
c----------------------------------------------------------------------
c
c     initialization 
      info = 0
      crr  = 0
c
c     call       
      crr = PRCRR(spot, strike, vol, rd, rf, tau, phi, ea, convention,
     &            quotation, greek, N) 
c
      RETURN
      END
