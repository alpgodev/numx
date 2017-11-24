c=======================================================================
c
c     subroutine MCFACT                                      
c
c     Monte-Carlo by factorial method
c
c-----------------------------------------------------------------------
      SUBROUTINE MCFACT ( USERPO, ndate, nasset, nfact, padata, pfdata,
     &                    H, N, T, nstep, rfr, indmu, mu,
     &                    iusdat, dusdat, iwork, dwork, 
     &                    iuwork, duwork, X, XF, price, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       USERPO : external subroutine (user payoff)             procedure
c       ndate  : number of dates                                 integer
c       nasset : number of assets                                integer
c       nfact  : number of factors                               integer
c       padata : assets prices, matrix (ndate*nasset)             double
c       pfdata : factors prices, matrix (ndate*nfact)             double
c       H      : horizon return in number of dates               integer
c       N      : number of simulations                           integer
c       T      : maturity (in years)                              double
c       nstep  : number of steps to maturity                     integer
c       rfr    : risk-free rate                                   double
c       indmu  : indicators, vector(nfact)                       integer
c                =1 : mu computed, =2 : mu gived
c       mu     : specified factor(s) drift, (nfact)               double
c       iusdat : user data for USERPO  vector                    integer
c                payoff function provided by the user
c       dusdat : user data for USERPO  vector                     double
c                payoff function provided by the user
c
c     WORKSPACE 
c       iwork  : ( 13*(nfact+1) )                                integer
c       dwork  : nfact*(2*ndate+3*nasset+7*nfact+2*nstep+43-H)
c                   + nasset*( 5*nasset+2*ndate+2*nstep+32-2*H )
c                   + ndate + 35  )                               double
c       iuwork : integer USERPO workspace                        integer
c       duwork : double USERPO workspace                          double
c
c     OUTPUT 
c       X[N]        : assets payoff (N)                           double
c       XF[nfact,N] : factors returns, (nfact,nsimul)             double
c       price       : payoff Monte-Carlo price                    double
c       info        : diagnostic argument                        integer
c
c     CALL   
c       YLMV    : copy a row of a vectorized matrix in a vector
c       HLRETM  : Computes the Horizon Log-RETurns on a Matrix
c       APTVOL  : Robusts alpha and beta (APT method) vectors of
c                 an universe for a given set of factors
c                 and volatility of regression error
c       CHOLMU  : computes the mean and the cholesky factor
c                 of the covariance matrix gamma of diffused prices
c                 derived from a log return distribution N(m,V)
c                 where m = vector such as mi=(mui-(sigmai)**2/2)*h
c                 and V=matrix h*gamma
c       PVX     : computing V*X = vector
c                 ( V vector(n), X scalar, gives V*X vector(n) )
c       MCFSIM  : computes one simulation of Monte-Carlo with factors
c       MV      : computing the mean of a vector
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL USERPO
c
c     arguments i/o
      INTEGER ndate, nasset, nfact, N, H, nstep, info
      INTEGER indmu(*), iusdat(*)
      DOUBLE PRECISION T, rfr, price
      DOUBLE PRECISION padata(*), pfdata(*), X(*), XF(*),dusdat(*),mu(*)
c
c     workspaces
      INTEGER iwork(*), iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER seed, nret, isim, ipo, imu, piw, pdw, pdra, pdrf, 
     &        pdalp, pdbet, pdvol, pdmu, pdcho, pdpif, pdpia
      DOUBLE PRECISION mpayoff, dt, sqrdt
#ifdef INTELFOR 
     
#else
      INTEGER time
#endif  
c
c     intrinsic function
      INTRINSIC DFLOAT, SQRT, EXP
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     dt and sqrt(dt)
      dt    = (T / nstep)
      sqrdt = SQRT(dt)
c
c     seed initialization for random
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
      pdrf  = pdra + (ndate-H)*nasset
c     pdrf  : pointer for factors returns matrix( (ndate-horiz)*nfact )
      pdalp = pdrf + (ndate-H)*nfact
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
c     Total size of dwork array = (nasset) + (nfact)
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
c     info = 0 (by construction)
      CALL YLMV ( ndate, nfact, pfdata, ndate, dwork(pdpif), info )
c     info = 0 (by construction)
c
c     log-returns (of assets)
      CALL HLRETM ( ndate, nasset, padata, H, dwork(pdra), info )
      IF (info .LT. 0) RETURN
c
c     log-returns (of factors)
      CALL HLRETM ( ndate, nfact, pfdata, H, dwork(pdrf), info )
      IF (info .LT. 0) RETURN
c
c     number of return(s)      
      nret = ndate - H
c
c     regression of assets returns on factors returns
      CALL APTVOL ( nret, nasset, nfact, dwork(pdra), dwork(pdrf),
     &              iwork(piw), dwork(pdw),
     &              dwork(pdalp), dwork(pdbet),
     &              dwork(pdvol), info )
      IF (info .LT. 0) RETURN
c
c     write alpha
c      open(unit=1,file='alpha.txt',status='unknown')
c      do i=1,nasset
c		write(1,*)  dwork(pdalp + i - 1)
c	end do   
c	close(unit=1)
c
c     write beta
c      open(unit=1,file='beta.txt',status='unknown')
c      do i=1,nfact*nasset
c		write(1,*)  dwork(pdbet + i - 1)
c	end do   
c	close(unit=1)
c
c     write volatility
c      open(unit=1,file='vol.txt',status='unknown')
c      do i=1,nasset
c		write(1,*)  dwork(pdvol + i - 1)
c	end do   
c	close(unit=1)
c
c     Choleski factorization and drifts
      CALL CHOLMU ( nret, nfact, H, dwork(pdrf),
     &              iwork(piw), dwork(pdw),
     &              dwork(pdmu), dwork(pdcho), info )
      IF (info .LT. 0) RETURN
c
c     expected return: mu 
      DO imu = 1,nfact
         IF( indmu(imu) .EQ. 2 ) THEN
            dwork(pdmu+imu-1) = mu(imu)
         ENDIF   
c          if ( indmu(imu).eq.2) then 
c            call RETMU ( givmu(imu), 100, pfdata(ndate*(imu-1)+1), 
c     &                   dwork(pdcho+(imu-1)*nfact + (imu-1)), matur,
c     &                   dwork(pdmu+imu-1), info)
c          endif
      ENDDO
c
c     mu*dt
      CALL PVX ( nfact, dwork(pdmu), dt, dwork(pdmu) )
c
c     write asset (initial)
c      open(unit=1,file='asset.txt',status='unknown')
c     do i=1,nasset
c		write(1,*)  dwork(pdpia + i - 1)
c	end do   
c	close(unit=1)
c
c     write factor (initial)
c      open(unit=1,file='factor.txt',status='unknown')
c      do i=1,nfact
c		write(1,*)  dwork(pdpif + i - 1)
c	end do   
c	close(unit=1)
c
c     write mu
c      open(unit=1,file='mu.txt',status='unknown')
c      do i=1,nfact
c		write(1,*)  dwork(pdmu + i - 1)
c	end do   
c	close(unit=1)
c
c     loop simulations
      ipo = 1
      DO isim = 1,N
        CALL MCFSIM ( USERPO, nasset, nfact,
     &                 dwork(pdalp), dwork(pdbet), dwork(pdvol),
     &                 dwork(pdpia), dwork(pdpif), H, T, nstep,
     &                 sqrdt, dwork(pdmu), dwork(pdcho),
     &                 iusdat, dusdat, dwork(pdw), iuwork, duwork,
     &                 X(isim), XF(ipo), info )
        IF (info .LT. 0) RETURN
        ipo = ipo + nfact
      ENDDO
c      
c     write X (payoff values)
c      open(unit=1,file='payoff.txt',status='unknown')
c      do i=1,N
c		write(1,*)  X(i)
c	end do   
c	close(unit=1)
c
c     write XF (factor(s) values)
c      open(unit=1,file='simuledFactor.txt',status='unknown')
c     do i=1,nfact*N
c		write(1,*) XF(i)
c	end do   
c	close(unit=1)
c
c     price = EXP(-r*T)*E[payoff]
      CALL MV ( N, X, mpayoff )
      price = EXP(-rfr*T) * mpayoff
c
      RETURN
      END
