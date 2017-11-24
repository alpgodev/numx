c=======================================================================
c
c     subroutine MCBASE                                      
c
c     Monte-Carlo method, basic version, input diffusion parmeters
c
c-----------------------------------------------------------------------
      SUBROUTINE MCBASE ( USERPO, p, mu, chol, prinit,
     &                    N, T, nstep, rfr,
     &                    iusdat, dusdat, dwork, iuwork,
     &                    duwork, X, price, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       USERPO : external payoff (provided by the user)        procedure
c       p      : dimension (p > 0)                               integer
c       mu     : drift, p-array                                   double
c       chol   : root square covariance (Choleski), matrix (p*p)  double
c       prinit : initial prices, p-array                          double
c       N      : number of simulations (max=100000)              integer
c       T      : maturity (in years)                              double
c       nstep  : number of steps to maturity  (nstep > 0)        integer
c       rfr    : risk-free rate                                   double
c       iusdat : user data for USERPO  vector                    integer
c                payoff function provided by the user
c       dusdat : user data for USERPO  vector                     double
c                payoff function provided by the user
c
c     WORKSPACE 
c       dwork  : ( p*(nstep + 5) )                                double
c       iuwork : integer workspace of payoff function            integer
c       duwork : double workspace of payoff function              double
c
c     OUTPUT 
c       X[N]   : payoff, N-array                                  double
c       price  : payoff Monte-Carlo price                         double
c       info   : diagnostic argument                             integer
c
c     CALL   
c       PVX    : computing V*X = vector ( V vector(n), X scalar, gives V*X vector(n) )
c       MCSIM  : one simulation of Monte-Carlo
c       USERPO : user payoff function
c       MV     : mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function
      EXTERNAL USERPO
c
c     i/o parameters
      INTEGER p, N, nstep, info
      INTEGER iusdat(*)
      DOUBLE PRECISION T, rfr, price
      DOUBLE PRECISION mu(*), prinit(*), X(*), chol(*), dusdat(*)
c
c     workspace
      INTEGER iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, i, pdw, pdtr, pdmud
      DOUBLE PRECISION dt, sqrdt, mpayoff
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c     intrinsic function
      INTRINSIC EXP
c
c-----------------------------------------------------------------------
c
c     initializations
      info  = 0        ! diagnostic argument
      dt    = T/nstep  ! time increment dt
      sqrdt = SQRT(dt) ! SQRT(dt)
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
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdmud = 1
c     pdmud : pointer for mu*deltat, ( p )
      pdtr  = pdmud + ( p )
c     pdtr  : pointer for trails matrix ((nstep + 1)*p),
      pdw   = pdtr + ( p*(nstep + 1) )
c     pdw   : pointer for MCSIM workspaces ( 3*p ) 
c
c     Total size of dwork array = (p) + (p*(nstep+1)) + (3*p)
c                               =  ( p*(nstep + 5) )
c
c-----------------------------------------------------------------------
c
c     mu*dt
      CALL PVX ( p, mu, dt, dwork(pdmud) )
c
c     loop simulation
      DO i = 1,N
c
c        Geometric Brownian Motion (Euler discretization)
c        S(t+1) = S(t)*[1 + mu*dt + chol(gamma)*dWt]
         CALL MCSIM ( p, nstep, sqrdt, prinit, dwork(pdmud), chol,
     &                dwork(pdw), dwork(pdtr) )
c
c        payoff evaluation
         CALL USERPO ( p, nstep+1, dwork(pdtr), T,
     &                 iusdat, dusdat, iuwork, duwork,
     &                 X(i), info )
         IF (info .LT. 0) RETURN
      ENDDO
c
c     mean price E[payoff]
      CALL MV ( N, X, mpayoff )
c
c     price = EXP(-r*T)*E[payoff]
      price = EXP( -rfr*T ) * mpayoff
c
      RETURN
      END
