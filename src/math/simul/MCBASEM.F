c=======================================================================
c
c     subroutine MCBASEM                                     
c
c     Monte-Carlo method, BASic version, input diffusion parmeters
c                         Multipayoff
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            USERPOM: entry point of an external subroutine
c                     multipayoff provided by the user         procedure
c            npotot : number of possible payoffs                 integer
c            indpo  : indicator of used payoffs, vectot(npotot)  integer
c                     used if indpo(i)=1
c            nasset : size of problem                            integer
c            mu     : drift, vector(nasset)                       double
c            chol   : choleski factor of covariance,
c                     matrix(nasset*nasset)                       double
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
c     WORKSPACE :
c            dwork  : vector( nasset*(nstep+5) )                  double
c            iuwork : user workspace for payoff function, vector integer
c            duwork : user workspace for payoff function, vector  double
c
c     OUTPUT :
c            mpay   : payoff, matrix(npotot,nsimul)               double
c            price  : price of each payoff, vectot(npotot)        double
c            info   : = 0 successful exit                        integer
c                     = -1 : nstep =< 0
c
c     CALL   :
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c        MCSIM   : computes one simulation of Monte-Carlo
c        USERPOM : Computes the multipayoff provided by the user
c        MLMI    : computing the mean of a row of a full matrix
c
c-----------------------------------------------------------------------
c
      subroutine MCBASEM ( USERPOM, npotot, indpo,
     &                     nasset, mu, chol, prinit,
     &                     nsimul, matur, nstep, rfrate,
     &                     liudat, iusdat, ldudat, dusdat,
     &                     dwork, iuwork, duwork,
     &                     mpay, price, info )
c
      implicit none
c
      external USERPOM
c
      integer npotot, nasset, nsimul, nstep, liudat, ldudat
      integer indpo(*), iusdat(liudat,*)
c
      double precision matur, rfrate
      double precision mu(*), prinit(*), mpay(npotot,*), price(*)
      double precision chol(*), dusdat(ldudat,*)
c
      integer iuwork(*)
      double precision dwork(*), duwork(*)
c
      integer info
c
      integer seed, isim, ipo
      integer pdw, pdtr, pdmud
      double precision deltat, sqdelt, moypay
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c     intrinsic function
      INTRINSIC SQRT, EXP
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     test if number of step > 0      
      IF (nstep .LE. 0) THEN
         info = -1
         RETURN
      ENDIF
c
      deltat = matur / nstep
      sqdelt = SQRT(deltat)
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
c     pdmud : pointer for mu*deltat, vector(nasset),
      pdtr  = pdmud + (nasset)
c     pdtr  : pointer for trails matrix((nstep+1)*nasset),
      pdw   = pdtr + (nasset*(nstep+1))
c     pdw   : pointer for local workspaces,
c             MCSIM needs (3*nasset)
c
c     Total size of dwork array = (nasset) + (nasset*(nstep+1))
c                               + (3*nasset)
c                               =  ( nasset *( (nstep+1) + 4 ) )
c
c             =  ( nasset*(nstep+5) )
c
c-----------------------------------------------------------------------
c
c     computing vector mu*deltat
      CALL PVX ( nasset, mu, deltat, dwork(pdmud) )
c
c     computing each simulation
      DO isim=1,nsimul
c
c        computing trails : S(t+1) = S(t) + dSt
c                           dSt = mu*dt + chol(gamma)*dWt
         CALL MCSIM ( nasset, nstep, sqdelt, prinit, dwork(pdmud), chol,
     &                dwork(pdw), dwork(pdtr) )
c
c        computing payoff
         call USERPOM ( npotot, indpo,
     &                  nasset, nstep+1, dwork(pdtr), matur,
     &                  liudat, iusdat, ldudat, dusdat, iuwork, duwork,
     &                  mpay(1,isim), info )
         IF (info .NE. 0) RETURN
      ENDDO
c
c     computing prices
      DO ipo = 1,npotot
         IF ( indpo(ipo) .EQ. 1 ) THEN
            CALL MLMI ( npotot, nsimul, mpay, ipo, moypay )
            price(ipo) = EXP(-rfrate*matur) * moypay
         ELSE
            price(ipo) = 0.0
         ENDIF
      ENDDO
c
      RETURN
      END
