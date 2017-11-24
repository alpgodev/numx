c=======================================================================
c
c     subroutine MCFACTM                                     
c
c     Monte-Carlo method with factors, multipayoff 
c
c-----------------------------------------------------------------------
      SUBROUTINE MCFACTM ( USERPOM, npotot, indpo, ndate, nasset, nfact,
     &                     padata, pfdata,
     &                     horiz, nsimul, matur, nstep, rfrate, epskat,
     &                     indmu, givmu,
     &                     liudat, iusdat, ldudat, dusdat,
     &                     iwork, dwork, iuwork, duwork,
     &                     mpay, pofact, price, info )
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
c            nfact  : number of factors                          integer
c            padata : prices assets during the estimation period
c                     matrix (ndate*nasset)                       double
c            pfdata : prics factors during the estimation period
c                     matrix (ndate*nfact)                        double
c            horiz  : horizon return in number of dates          integer
c            nsimul : number of simulations                      integer
c            matur  : maturity (in years)                         double
c            nstep  : number of steps to maturity                integer
c            rfrate : risk-free rate                              double
c            epskat : Kato sensibility factor                     double
c            indmu  : indicators for mu, vector(nfact)           integer
c                     =1 : mu computed, =2 : mu gived
c            givmu  : gived mu, vector(nfact)                     double
c            liudat : size of integer user data for USERPOM      integer
c            iusdat : user data for USERPOM,  matrix(liudat,*)   integer
c                     payoff function provided by the user
c            ldudat : size of double prec. user data for USERPOM integer
c            dusdat : user data for USERPOM, matrix(ldudat,*)     double
c                     payoff function provided by the user
c
c     WORKSPACE :
c            iwork  : vector( 13*(nfact+1) )                     integer
c            dwork  : vector(  
c                     nfact*(2*ndate+3*nasset+7*nfact+2*nstep+43-horiz)
c                   + nasset*( 5*nasset+2*ndate+2*nstep+32-2*horiz )
c                   + ndate + 35  )                               double
c            iuwork : user workspace for payoff function, vector integer
c            duwork : user workspace for payoff function, vector  double
c
c     OUTPUT :
c            mpay   : payoff, matrix(npotot,nsimul)               double
c            pofact : factors payoff, matrix(nfact,nsimul)        double
c            price  : price of each payoff, vectot(npotot)        double
c            info   : = 0 successful exit                        integer
c                     = -1 : nstep =< 0
c
c     CALL   :
c        YLMV    : copy a row of a vectorized matrix in a vector
c        HLRETM  : Computes the Horizon Log-RETurns on a Matrix
c        APTVOL  : Robusts alpha and beta (APT method) vectors of
c                  an universe for a given set of factors
c                  and volatility of regression error
c        CHOLMU  : computes the mean and the cholesky factor
c                  of the covariance matrix gamma of diffused prices
c                  derived from a log return distribution N(m,V)
c                  where m = vector such as mi=(mui-(sigmai)**2/2)*h
c                  and V=matrix h*gamma
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c        MCFSIMM : computes one simulation of Monte-Carlo with factors,
c                  multipayoff 
c        MLMI    : computing the mean of a row of a full matrix
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL USERPOM
c
c     arguments i/o
      INTEGER npotot, ndate, nasset, nfact, nsimul, horiz, nstep,
     &        liudat, ldudat, info
      INTEGER indpo(*), indmu(*), iusdat(*)
      DOUBLE PRECISION matur, rfrate, epskat
      double precision padata(*), pfdata(*), mpay(npotot,*), price(*),
     &                 pofact(*), dusdat(*), givmu(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, nret, isim, ipo, ifact, imu, piw, pdw, pdra, 
     &        pdrf, pdalp, pdbet, pdvol, pdmu, pdcho, pdpif, pdpia
      DOUBLE PRECISION moypay, deltat, sqdelt
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
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
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for local workspaces,
c             APTVOL needs 13*(nfact+1)
c             CHOLMU needs (12*nfact)
c             so, union of spaces : 13*(nfact+1)
c
c     Total size of iwork array = 13*(nfact+1)
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdpia = 1
c     pdpia : pointer for initial prices of assets, vector(nasset)
      pdpif = pdpia + (nasset)
c     pdpif : pointer for initial prices of factors, vector(nfact)
      pdra  = pdpif + (nfact)
c     pdra  : pointer for assets returns matrix( (ndate-horiz)*nasset )
      pdrf  = pdra + (ndate-horiz)*nasset
c     pdrf  : pointer for factors returns matrix( (ndate-horiz)*nfact )
      pdalp = pdrf + (ndate-horiz)*nfact
c     pdalp : pointer for alphas, vector(nasset)
      pdbet = pdalp + (nasset)
c     pdbet : pointer for betas, matrix(nfact*nasset)
      pdvol = pdbet + (nfact*nasset)
c     pdvol : pointer for volatilities, vector(nasset)
      pdmu  = pdvol + (nasset)
c     pdmu  : pointer for drifts mu of factors, vector(nfact)
      pdcho = pdmu + (nfact)
c     pdcho : pointer for Choleski factor, matrix(nfact*nfact)
      pdw   = pdcho + (nfact*nfact)
c     pdw   : pointer for local workspaces,
c             APTVOL needs  nfact*( ndate + 2*nasset + 6*nfact + 41 )
c                            + ( ndate*nasset + ndate + 2*nasset + 35 )
c             CHOLMU needs  nasset*( 5*nasset+nret+29 )
c                           = nasset*( 5*nasset+ndate-horiz+29 )
c             MCFSIM needs  (2*nstep+4)*nfact) + ((2*nstep+1)*nasset)
c             so, union of spaces :
c                      nfact*( ndate + 2*nasset + 6*nfact + 41 )
c                      + ( ndate*nasset + ndate + 2*nasset + 35 )
c                      +  nasset*( 5*nasset - horiz + 27 )
c                      + (2*nstep*nfact) + (2*nstep*nasset) 
c
c      Total size of dwork array = (nasset) + (nfact)
c                      + (ndate-horiz)*nasset + (ndate-horiz)*nfact
c                      + (nasset) + (nfact*nasset) + (nasset) + (nfact)
c                      + (nfact*nfact)
c                      + nfact*( ndate + 2*nasset + 6*nfact + 41 )
c                      + ( ndate*nasset + ndate + 2*nasset + 35 )
c                      +  nasset*( 5*nasset - horiz + 27 )
c                      + ( nstep*nfact ) 
c
c             = nfact*( 1+ndate-horiz+nasset+1+nfact
c                       +ndate+2*nasset+6*nfact+41+2*nstep )
c             + nasset*( 1+ndate-horiz+1+1+ndate+2
c                        +(5*nasset+2*nstep-horiz+27) )
c             + ( ndate + 35 )
c
c             = nfact*( 2*ndate +3*nasset +7*nfact +2*nstep +43 -horiz )
c             + nasset*( 5*nasset + 2*ndate + 2*nstep + 32 - 2*horiz )
c             + ( ndate + 35 )
c
c-----------------------------------------------------------------------
c
c     saves initial prices of assets and factors (last price of history)
      CALL YLMV ( ndate, nasset, padata, ndate, dwork(pdpia), info )
      CALL YLMV ( ndate, nfact, pfdata, ndate, dwork(pdpif), info )
c     info = 0 by construction
c
c     Log-returns (assets)
      CALL HLRETM ( ndate, nasset, padata, horiz, dwork(pdra), info)
      IF (info .LT. 0) RETURN
      nret = ndate - horiz
c
c     Log-returns (factors)
      CALL HLRETM ( ndate, nfact, pfdata, horiz,
     &              dwork(pdrf), info )
      IF (info .LT. 0) RETURN
c
c     regression of assets on factor(s)
      CALL APTVOL ( nret, nasset, nfact, dwork(pdra), dwork(pdrf),
     &              epskat, iwork(piw), dwork(pdw),
     &              dwork(pdalp), dwork(pdbet),
     &              dwork(pdvol), info )
      IF (info .LT. 0) RETURN
c
c     Choleski factorization and drifts
      CALL CHOLMU ( nret, nfact, horiz, dwork(pdrf),
     &              iwork(piw), dwork(pdw),
     &              dwork(pdmu), dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
c     if index = 2 then mu = specified by user
      DO imu = 1,nfact
         IF (indmu(imu) .EQ. 2) dwork(pdmu+imu-1) = givmu(imu)
      ENDDO
c
c     mu*deltat
      CALL PVX ( nfact, dwork(pdmu), deltat, dwork(pdmu) )
c
c     loop (simulation)
      ifact = 1
      DO isim = 1,nsimul
c
         CALL MCFSIMM ( USERPOM, npotot, indpo, nasset, nfact,
     &                  dwork(pdalp), dwork(pdbet), dwork(pdvol),
     &                  dwork(pdpia), dwork(pdpif), horiz, matur, nstep,
     &                  sqdelt, dwork(pdmu), dwork(pdcho),
     &                  liudat, iusdat, ldudat, dusdat,
     &                  dwork(pdw), iuwork, duwork,
     &                  mpay(1,isim), pofact(ifact), info )
         IF (info .LT. 0) RETURN
c
         ifact = ifact + nfact
      ENDDO
c
c     prices
      DO ipo=1,npotot
         IF (indpo(ipo) .EQ. 1) THEN
            CALL MLMI ( npotot, nsimul, mpay, ipo, moypay )
            price(ipo) = EXP(-rfrate*matur) * moypay
         ELSE
            price(ipo) = 0.0
         ENDIF
      ENDDO
c
      RETURN
      END
