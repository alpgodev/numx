c=======================================================================
c
c     subroutine MCBCORM                                     
c
c     Monte-Carlo method, Basic version, input correlation matrix
c                         Multipayoff
c
c-----------------------------------------------------------------------
      SUBROUTINE MCBCORM ( USERPOM, npotot, indpo,
     &                     nasset, matcor, mu, volat, prinit,
     &                     nsimul, matur, nstep, rfrate,
     &                     liudat, iusdat, ldudat, dusdat,
     &                     iwork, dwork, iuwork, duwork,
     &                     mpay, price, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            USERPOM: entry point of an external subroutine
c                     multipayoff provided by the user         procedure
c            npotot : number of possible payoffs                 integer
c            indpo  : indicator of used payoffs, vectot(npotot)  integer
c                     used if indpo(i)=1
c            nasset : size of problem                            integer
c            matcor : correlation, matrix(nasset*nasset)          double
c            mu     : drift, vector(nasset)                       double
c            volat  : volatilities, vector(nasset)                double
c            prinit : initial prices, vector(nasset)              double
c            nsimul : number of simulations                      integer
c            matur  : maturity (in years)                         double
c            nstep  : number of steps to maturity                integer
c            rfrate : risk-free rate                              double
c            liudat : size of integer user data for USERPOM      integer
c            iusdat : user data for USERPOM,  matrix(liudat,*)   integer
c                     payoff function provided by the user
c            ldudat : size of double prec. user data for USERPOM integer
c            dusdat : user data for USERPOM, matrix(ldudat,*)     double
c                     payoff function provided by the user
c
c     WORKSPACE 
c            iwork  : vector( 12*nasset )                        integer
c            dwork  : vector( nasset*(6*nasset+nstep+27) )        double
c            iuwork : user workspace for payoff function, vector integer
c            duwork : user workspace for payoff function, vector  double
c
c     OUTPUT 
c            mpay   : payoff, matrix(npotot,nsimul)               double
c            price  : price of each payoff, vectot(npotot)        double
c            info   : = 0 successful exit                        integer
c                     = -1 : nstep =< 0
c
c     CALL   
c        CORCOV  : conversion correlation/covariance matrix
c        RCHO    : Robust Choleski factorization
c        MCBASEM : Monte-Carlo method, BASic version, Multiple payoff,
c                  input diffusion parmeters
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function
      EXTERNAL USERPOM
c
c     arguments i/o
      INTEGER info, npotot, nasset, nsimul, nstep, liudat, ldudat
      INTEGER indpo(*), iusdat(liudat,*)
      DOUBLE PRECISION matur, rfrate
      DOUBLE PRECISION mu(*), volat(*), prinit(*), matcor(*), 
     &                 dusdat(ldudat,*), mpay(npotot,*), price(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, piw, pdw, pdcho, pdcov
      DOUBLE PRECISION deltat
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c     intrinsic function
      INTRINSIC DFLOAT
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      IF (nstep .LE. 0) THEN
         info = -1
         RETURN
      ENDIF
      deltat = (matur / nstep)
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
c     piw   : pointer for local workspaces,
c             RCHO   needs (12*nasset)
c
c     Total size of iwork array = (12*nasset)
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdcov = 1
c     pdcov : pointer for covariance matrix(nasset*nasset),
c             so (nasset*nasset) more,
      pdcho = pdcov + (nasset*nasset)
c     pdcho : pointer for choleski matrix(nasset*nasset),
c             so (nasset*nasset) more,
      pdw   = pdcho + (nasset*nasset)
c     pdw   : pointer for local workspaces,
c             RCHO    needs ( nasset*(4*nasset+27) )
c             MCBASE needs ( nasset*(nstep+5) )
c             so ( nasset*(4*nasset+nstep+27) ) (union of workspaces)
c
c     Total size of dwork array = (nasset*nasset) + (nasset*nasset)
c                                    + ( nasset*(4*nasset+nstep+27) )
c
c             = ( nasset *( 6*nasset + nstep + 27 ) )
c
c-----------------------------------------------------------------------
c
c     transformation correlation matrix in covariance matrix
      CALL CORCOV ( nasset, matcor, volat, dwork(pdcov) )
c
c     Choleski factorization
      CALL RCHO ( nasset, dwork(pdcov), iwork(piw), dwork(pdw),
     &            dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
c     computing Monte-Carlo
      CALL MCBASEM ( USERPOM, npotot, indpo,
     &               nasset, mu, dwork(pdcho), prinit,
     &               nsimul, matur, nstep, rfrate,
     &               liudat, iusdat, ldudat, dusdat,
     &               dwork(pdw), iuwork, duwork,
     &               mpay, price, info )
c
      RETURN
      END
