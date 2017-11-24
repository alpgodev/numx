c=======================================================================
c
c     subroutine MCBHISM                                     
c
c     Monte-Carlo method, basic version,
c     input assets prices history, Multipayoff
c
c-----------------------------------------------------------------------
      SUBROUTINE MCBHISM ( USERPOM, npotot, indpo,
     &                     ndate, nasset, padata, horiz,
     &                     nsimul, matur, nstep, rfrate, indmu, givmu,
     &                     liudat, iusdat, ldudat, dusdat,
     &                     iwork, dwork, iuwork, duwork,
     &                     mpay, price, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            USERPOM: entry point of an external subroutine
c                     multipayoff provided by the user         procedure
c            npotot : number of possible payoffs                 integer
c            indpo  : indicator of used payoffs, vectot(npotot)  integer
c                     used if indpo(i)=1
c            ndate  : number of dates                            integer
c            nasset : number of assets                           integer
c            padata : prices assets during the estimation period
c                     matrix (ndate*nasset)                       double
c            horiz  : horizon return in number of dates          integer
c            nsimul : number of simulations                      integer
c            matur  : maturity (in years)                         double
c            nstep  : number of steps to maturity                integer
c            rfrate : risk-free rate                              double
c            indmu  : indicators for mu, vector(nasset)          integer
c                     =1 : mu computed, =2 : mu gived
c            givmu  : gived mu, vector(nasset)                    double
c            liudat : size of integer user data for USERPOM      integer
c            iusdat : user data for USERPOM,  matrix(liudat,*)   integer
c                     payoff function provided by the user
c            ldudat : size of double prec. user data for USERPOM integer
c            dusdat : user data for USERPOM, matrix(ldudat,*)     double
c                     payoff function provided by the user
c
c     WORKSPACE :
c            iwork  : vector( 12*nasset )                        integer
c            dwork  : vector(
c                     nasset*( 5*nasset+2*ndate+2*nstep-2*horiz+31 )
c                           )                                     double
c            iuwork : user workspace for payoff function, vector integer
c            duwork : user workspace for payoff function, vector  double
c
c     OUTPUT :
c            mpay   : payoff, matrix(npotot,nsimul)               double
c            price  : price of each payoff, vectot(npotot)        double
c            cinfo  : path of error                         character*80
c            info   : = 0 successful exit                        integer
c                     = -1 : nstep =< 0
c                     = else : see MCBASE guide
c
c     CALL
c        YLMV    : copy a row of a vectorized matrix in a vector
c        HLRETM  : Computes the Horizon Log-RETurns on a Matrix
c        CHOLMU  : computes the mean and the cholesky factor
c                  of the covariance matrix gamma of diffused prices
c                  derived from a log return distribution N(m,V)
c                  where m = vector such as mi=(mui-(sigmai)**2/2)*h
c                  and V=matrix h*gamma
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
      INTEGER info, npotot, ndate, nasset, nsimul, nstep, horiz,
     &        liudat, ldudat
      INTEGER indpo(*), indmu(*), iusdat(liudat,*)
      DOUBLE PRECISION rfrate, matur
      DOUBLE PRECISION padata(*), dusdat(ldudat,*), givmu(*), 
     &                 mpay(npotot,*), price(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
      INTEGER seed, nret, imu, piw, pdw, pdpini, pdra, pdmu, pdcho
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
c             CHOLMU needs (12*nfact)
c
c     Total size of iwork array = (12*nasset)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdpini= 1
c     pdpini: pointer for initial prices of assets, vector(nasset)
      pdra  = pdpini + (nasset)
c     pdra  : pointer for assets returns matrix( (ndate-horiz)*nasset )
      pdmu  = pdra + (ndate-horiz)*nasset
c     pdmu  : pointer for drifts mu of factors, vector(nasset)
      pdcho = pdmu + (nasset)
c     pdcho : pointer for Choleski factor, matrix(nasset*nasset)
      pdw   = pdcho + (nasset*nasset)
c     pdw   : pointer for local workspaces,
c             CHOLMU needs  nasset*( 5*nasset+nret+29 )
c                           = nasset*( 5*nasset+ndate-horiz+29 )
c             MCBASE needs ( nasset*(nstep+5) )
c             so (union of workspaces):
c             nasset*( 5*nasset+ndate-horiz+nstep+29 )
c
c     Total size of dwork array = (nasset) + (ndate-horiz)*nasset
c                               + (nasset)
c                               + (nasset*(nstep+1))
c                               +  nasset*( 5*nasset+ndate-horiz+nstep+29 )
c             = nasset*( 1 + ndate-horiz + 1 + nstep+1
c                        + 5*nasset+ndate-horiz+nstep+29  )
c
c             = ( nasset*( 5*nasset+2*ndate+2*nstep-2*horiz+31 ) )
c
c-----------------------------------------------------------------------
c
c     saves initial prices of assets (last price of history)
      CALL YLMV ( ndate, nasset, padata, ndate, dwork(pdpini), info )
c     info = 0 by construction
c
c     computes Log-returns (assets)
      CALL HLRETM ( ndate, nasset, padata, horiz,
     &              dwork(pdra), info )
      IF (info .LT. 0) RETURN
      nret = ndate - horiz
c
c     Choleski factorization (volatility) and drifts
      CALL CHOLMU ( nret, nasset, horiz, dwork(pdra),
     &              iwork(piw), dwork(pdw),
     &              dwork(pdmu), dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
      DO imu = 1,nasset
         IF ( indmu(imu).eq.2 ) dwork(pdmu+imu-1) = givmu(imu)
      ENDDO
c
c     computing Monte-Carlo
      CALL MCBASEM ( USERPOM, npotot, indpo,
     &               nasset, dwork(pdmu), dwork(pdcho),
     &               dwork(pdpini), nsimul, matur, nstep, rfrate,
     &               liudat, iusdat, ldudat, dusdat,
     &               dwork(pdw), iuwork, duwork,
     &               mpay, price, info )
      RETURN
      END
