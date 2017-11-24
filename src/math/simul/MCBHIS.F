c=======================================================================
c
c     subroutine MCBHIS                                     
c
c     Monte-Carlo method, basic version, input assets prices HIStory
c
c-----------------------------------------------------------------------
      SUBROUTINE MCBHIS ( USERPO, ndate, P, padata, H,
     &                    N, T, nstep, rfr, indmu, mu,
     &                    liudat, iusdat, ldudat, dusdat, iwork, dwork, 
     &                    iuwork, duwork, X, price, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       USERPO   : external subroutine (user payoff)          procedure
c       ndate    : number of dates                              integer
c       P        : number of asset(s)                           integer
c       padata   : assets(s) prices (ndate*P)                    double
c       H        : horizon return in number of dates            integer
c       N        : number of simulations                        integer
c       T        : maturity (in years)                           double
c       nstep    : number of steps to maturity                  integer
c       rfr      : risk-free rate                                double
c       indmu[P] : indicators for mu (P)                        integer
c                  =1 : mu computed, =2 : mu gived
c       mu[P]    : specified mu (P)                              double
c       liudat   : size of the integer workspace                integer
c       iusdat   : integer user data for USERPO function        integer
c       ldudat   : size of the integer workspace                integer
c       dusdat   : double user data for USERPO function          double
c
c     WORKSPACE 
c       iwork    : ( 12*P )                                     integer
c       dwork    : P*( 5*P + 2*ndate + 2*nstep - 2*H + 31 )      double
c       iuwork   : integer USERPO workspace                     integer
c       duwork   : double USERPO workspace                       double
c
c     OUTPUT 
c       X[N]     : payoff (N)                                    double
c       price    : payoff Monte-Carlo price                      double
c       info     : diagnostic argument                          integer
c
c     CALL   
c       YLMV    : copy a row of a vectorized matrix in a vector
c       HLRETM  : Computes the Horizon Log-RETurns on a Matrix
c       CHOLMU  : computes the mean and the cholesky factor
c                 of the covariance matrix gamma of diffused prices
c                 derived from a log return distribution N(m,V)
c                 where m = vector such as mi=(mui-(sigmai)**2/2)*h
c                 and V=matrix h*gamma
c       MCBASE  : Monte-Carlo method, BASic version
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function (payoff)
      EXTERNAL USERPO
c
c     i/o parameters
      INTEGER ndate, P, N, nstep, H, info, indmu(*), iusdat(*)
      INTEGER liudat, ldudat
      DOUBLE PRECISION rfr, T, price
      DOUBLE PRECISION padata(*), X(*), dusdat(*), mu(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, nret, imu, piw, pdw, pdpini, pdra, pdmu, pdcho
      DOUBLE PRECISION dt
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     time increment
      dt = T / nstep
c
c     seed for random
#ifdef INTELFOR 
      CALL TIME(seed)
      CALL RANDOM_SEED(seed)
#else
      seed = time()
      CALL srand(seed)
#endif      
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for CHOLMU workspaces (12*P)
c
c     Total size of iwork array = ( 12*P )
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdpini= 1
c     pdpini: pointer for initial prices of assets, vector(P)
      pdra  = pdpini + ( P )
c     pdra  : pointer for assets returns matrix( (ndate-H)*P )
      pdmu  = pdra + ((ndate-H)*P)
c     pdmu  : pointer for drifts mu of factors, vector(P)
      pdcho = pdmu + ( P )
c     pdcho : pointer for Choleski factor, matrix(P*P)
      pdw   = pdcho + ( P*P )
c     pdw   : pointer for local workspaces,
c             CHOLMU needs    P*( 5*P + nret + 29 )
c                           = P*( 5*P + ndate - H + 29 )
c             MCBASE needs ( P*(nstep+5) )
c             so (union of workspaces):
c             P*( 5*P+ndate - H + nstep + 29 )
c
c     Total size of dwork array = (P) + (ndate - H)*P
c                            + (P)
c                            + (P*(nstep + 1))
c                            +  P*( 5*P + ndate - H + nstep + 29 )
c             = ( P*( 5*P + 2*ndate + 2*nstep - 2*H + 31 ) )
c
c-----------------------------------------------------------------------
c
c     saves initial prices of assets (last price of history)
      CALL YLMV ( ndate, P, padata, ndate, dwork(pdpini), info )
c     info = 0 by construction
c
c     log-returns (of assets)
      CALL HLRETM ( ndate, P, padata, H, dwork(pdra), info )
      IF ( info .LT. 0 ) RETURN
      nret = ndate - H
c
c     Choleski factorization (volatility) and drifts
      CALL CHOLMU ( nret, P, H, dwork(pdra),
     &              iwork(piw), dwork(pdw),
     &              dwork(pdmu), dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
      DO imu = 1,P
         IF ( indmu(imu) .EQ. 2 ) dwork(pdmu+imu-1) = mu(imu)
      ENDDO
c
c     Monte-Carlo (basic form)
      CALL MCBASE ( USERPO, P, dwork(pdmu), dwork(pdcho),
     &              dwork(pdpini), N, T, nstep, rfr,
     &              liudat, iusdat, ldudat, dusdat, dwork(pdw), iuwork, 
     &              duwork, X, price, info )
      RETURN
      END
