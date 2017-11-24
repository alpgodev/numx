c=======================================================================
c
c     subroutine MCBCOR                                   
c
c     Monte-Carlo method, Basic version, input correlation matrix
c
c-----------------------------------------------------------------------
      SUBROUTINE MCBCOR ( USERPO, P, corr, mu, vol, prinit, N, T, 
     &                    nstep, rfr,
     &                    liudat, iusdat, ldudat, dusdat, iwork, dwork, 
     &                    iuwork, duwork, X, price, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       USERPO    : external subroutine (user payoff)          procedure
c       P         : number of asset(s)                           integer
c       corr[P,P] : correlation matrix (P*P)                      double
c       mu[P]     : diffusion drift (P)                           double
c       vol[P]    : volatilities (P)                              double
c       prinit[P] : initial prices (P)                            double
c       N         : number of simulations                        integer
c       T         : maturity (in years)                           double
c       nstep     : number of steps to maturity                  integer
c       rfr       : risk-free rate                                double
c       liudat   : size of the integer workspace                 integer
c       iusdat : user data for USERPO  vector                    integer
c                payoff function provided by the user
c       ldudat   : size of the integer workspace                 integer
c       dusdat : user data for USERPO  vector                     double
c                payoff function provided by the user
c
c     WORKSPACE 
c       iwork     : ( 12*P )                                     integer
c       dwork     : ( P*(6*P + nstep + 27) )                      double
c       iuwork    : integer USERPO workspace                     integer
c       duwork    : double USERPO workspace                       double
c
c     OUTPUT 
c       X[N]      : payoff (N)                                    double
c       price     : payoff Monte-Carlo price                      double
c       info      : diagnostic argument                          integer
c
c     CALL   
c       CORCOV    : correlation matrix -> covariance matrix
c       RCHO      : robust Choleski factorization
c       MCBASE    : Monte-Carlo method, BASic version
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function (payoff)
      EXTERNAL USERPO
c
c     i/o parameters
      INTEGER P, N, nstep, info, iusdat(*), liudat, ldudat 
      DOUBLE PRECISION T, rfr, price
      DOUBLE PRECISION mu(*), vol(*), prinit(*), X(*), corr(*),
     & dusdat(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, piw, pdw, pdcho, pdcov
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c-----------------------------------------------------------------------
c
c     initializations
      info   = 0
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
c     piw   : pointer for RCHO workspaces, (12*P)
c
c     Total size of iwork array = ( 12*P )
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdcov = 1
c     pdcov : pointer for covariance matrix (P*P),
      pdcho = pdcov + ( P*P )
c     pdcho : pointer for choleski matrix(P*P),
      pdw   = pdcho + ( P*P )
c     pdw   : pointer for local workspaces,
c             RCHO    needs ( P*(4*P + 27) )
c             MCBASE  needs ( P*(nstep + 5) )
c             so ( P*(4*P + nstep + 27) ) (union of workspaces)
c
c     Total size of dwork array = ( P*( 6*P + nstep + 27 ) )
c
c-----------------------------------------------------------------------
c
c     correlation matrix -> covariance matrix
      CALL CORCOV ( P, corr, vol, dwork(pdcov) )
c
c     Choleski factorization
      CALL RCHO ( P, dwork(pdcov), iwork(piw), dwork(pdw),
     &            dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
c     Monte-Carlo (basic form)
      CALL MCBASE ( USERPO, P, mu, dwork(pdcho), prinit,
     &              N, T, nstep, rfr,
     &              liudat, iusdat, ldudat, dusdat, dwork(pdw), iuwork, 
     &              duwork, X, price, info )
c
      RETURN
      END
