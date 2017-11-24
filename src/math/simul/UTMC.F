c=======================================================================
c
c     Utility Functions                                    
c
c     subroutines for Monte-Carlo method
c
c-----------------------------------------------------------------------
c
c        GBMSTEP : one step of a p-dim Geometric Brownian Motion (BGM)
c        CHOLMU  : computes the mean and the cholesky factor
c                  of the covariance matrix gamma of diffused prices
c                  derived from a log return distribution N(m,V)
c                  where m = vector such as mi=(mui-(sigmai)**2/2)*h
c                  and V=matrix h*gamma
c        MCSIM   : one simulation of Monte-Carlo
c        MCFSIM  : one simulation of Monte-Carlo with factors
c        MCFSIMM : one simulation of Monte-Carlo with factors, multipayoff
c        MURET   : diffusion parameter mu from return to maturity of the asset
c                  signature to be used by DICHOT
c                  (dichotomy of a monotonous function)
c
c-----------------------------------------------------------------------
c
c     subroutine GBMSTEP                                    
c
c     Generates one step of a multidimensional Geometric Brownian Motion
c     by Euler discretization method
c
c     dS(t) = S(t)*[mu*dt + sigma*dW(t)]
c
c-----------------------------------------------------------------------
      SUBROUTINE GBMSTEP ( n, mudt, chol, sqdt, dwork, step )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n         : number of elements                           integer
c       mudt[n]   : mu(i)*deltat, (n)                             double
c       chol[n,n] : choleski factor of covariance (n*n)           double
c       sqdt      : sqrt of dt (discret. step in years)           double
c
c     WORKSPACE
c       dwork     : ( 2*n )                                       double
c
c     OUTPUT 
c       step[n]   : generated step (n)                            double
c
c     CALL   
c       PMLV      : M*V = vector (M low triangular matrix(n*n), V vector(n),
c                   gives M*V vector(n))
c       SVVX      : V1 + V2*X, V1 and V2 vectors(n), X scalar, give vector(n)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n
      DOUBLE PRECISION sqdt
      DOUBLE PRECISION mudt(*), chol(*), step(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdb, pdc
      DOUBLE PRECISION mean, std
      PARAMETER ( mean = 0.0d0, std = 1.0d0 )
c
c-----------------------------------------------------------------------
c
c     initialization
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdb   = 1
c     pdb   : pointer for noise vector(nelem), so (n) more,
      pdc   = pdb + ( n )
c     pdc   : pointer for work vector(nelem), so (n) more,
c
c     Total size of dwork array = (n) + (n) = ( 2*n )
c
c-----------------------------------------------------------------------
c
c     generate n-Gaussian number N(m,std)
      CALL IVXRN ( n, mean, std, dwork(pdb) )
c
c     sigma*dW, vector(n)
      CALL PMLV ( n, chol, dwork(pdb), dwork(pdc) ) 
c
c     [mu*dt] + [sigma*dW*sqrt(t)], vector(n)
      CALL SVVX (  n, mudt, dwork(pdc), sqdt, step ) 
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine CHOLMU                                      
c
c     This function computes the mean and the cholesky factor of the cov. 
c     matrix of diffused prices derived from a return distribution N(m,V)
c     where m = vector such as mi=(mui-0.5*(sigmai)**2)*h
c     and V = matrix h*gamma
c
c-----------------------------------------------------------------------
      SUBROUTINE CHOLMU ( n, p, h, X, iwork, dwork, mu, chol, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n         : number of returns (n > 1)                    integer
c       p         : number of asset(s) (                         integer
c       h         : horizon return in number of dates            integer
c       X[n,p]    : returns, matrix (n*p)                         double
c
c     WORK
c       iwork     : ( 12*p )                                     integer
c       dwork     : p*( 5*p + n + 29 )                            double
c
c     OUTPUT 
c       mu[p]     : drifts of the process(es) (p)                 double
c       chol[p,p] : Choleski, square root covariance (p*p)        double
c       info      : diagnostic argument                          integer
c
c     CALL   
c       COVM      : covariance matrix and mean vector (empirical)
c       PMX       : M*X = matrix ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c       CMDV      : converting the diagonal of a matrix(n*n) in vector(n)
c       SVXVY     : computing V1*X + V2*Y 
c                   V1 and V2 vectors(n), X and Y scalars, give vector(n)
c       RCHO      : robust Choleski factorization
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, h, info
      DOUBLE PRECISION X(*), mu(*), chol(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdm, pdcov, pdvar, pdw, piw, nbdy
      DOUBLE PRECISION invh, demi
      PARAMETER ( demi = 0.5, nbdy = 252 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     test if horizon >= 1      
      IF (h .LE. 0) THEN
         info = -401
         RETURN
      ENDIF
c
c     nb. days in year / horizon      
      invh = float(nbdy)/h
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piw = 1
c     piw   : pointer for RCHO workspace who needs (12*nmat)
c
c     Total size of iwork array : (12*nmat)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdm   = 1
c     pdm   : pointer for means, vector(p)
      pdvar = pdm + ( p )
c     pdvar : pointer for diagonalof cov, vector(p)
      pdcov = pdvar + ( p )
c     pdcov : pointer for covariance matrix, matrix(p*p)
      pdw   = pdcov + ( p*p )
c     pdw   : pointer for workspaces,
c             COVMVM needs ( n*p + p)
c             RCHO needs p*( 4*p + 27 )
c             so p*( 5*p + n + 27 ) more
c
c     Total size of dwork array = p*( 5*p + n + 29 )
c
c-----------------------------------------------------------------------
c
c     mean and covariance matrix of log-returns 
      CALL COVMVM ( n, p, X, dwork(pdw), dwork(pdm), dwork(pdcov), 
     &            info)
c
c     write covariance matrix
c	call IMPRMC ( 10, 'cov. matrix:', p, dwork(pdcov) )
c
c     divides covariance matrix by horizon
      CALL PMX ( p, p, dwork(pdcov), invh, dwork(pdcov) )
c
c     extract variances of covariance matrix (diagonal)
      CALL CMDV ( p, dwork(pdcov), dwork(pdvar) )
c
c     mu = mean/horizon + variance/2.
      CALL SVXVY ( p, dwork(pdm), invh, dwork(pdvar), demi, mu )
c
c     Choleski factor
      CALL RCHO ( p, dwork(pdcov), iwork(piw), dwork(pdw), 
     &            chol, info )
c
c     write Cholesky matrix
c     open(unit=1,file='nelem.txt',status='unknown')
c	write(1,*)  p
c	close(unit=1)
c	call IMPRMC ( 10, 'Chol. matrix:', p, chol )     
c     
      RETURN
      END
c
c=======================================================================
c
c     subroutine  MCSIM                                     
c
c     Computes one simulation of Monte-Carlo
c
c-----------------------------------------------------------------------
      SUBROUTINE MCSIM ( p, N, sqrdt, pinit, mudt, chol, dwork, trails )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p         : number of element(s) (p > 0)                 integer
c       N         : number of steps to maturity                  integer
c       sqrdt     : sqrt of deltat (discret. step in years)       double
c       pinit[p]  : initial values, (p)                           double
c       mudt[p]   : mu(i)*deltat, (p)                             double
c       chol[p,p] : choleski factor of covariance (p*p)           double
c
c     WORKSPACE 
c       dwork     : ( 3*p )                                       double
c
c     OUTPUT 
c       trails[N+1,p] : trails, matrix((N+1)*p)                   double
c
c     CALL   
c        YVLM     : copy a vector in a row of a vectorized matrix
c        GBMSTEP  : multidimensional GBM: dSt = mu*dt + chol(gamma)*dWt
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, N, info
      DOUBLE PRECISION sqrdt
      DOUBLE PRECISION pinit(*), mudt(*), chol(*), trails((N+1),*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdstep, pdw, i, k, ks
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdstep= 1
c     pdstep: pointer for steps on each factor, vector(p)
      pdw   = pdstep + (p)
c     pdw   : pointer for GBMSTEP workspaces (2*p)
c
c     Total size of dwork array = 3*p
c
c-----------------------------------------------------------------------
c
c     initialization of factors trails at initial price
      CALL YVLM ( N+1, p, trails, pinit, 1, info )
c
c     trails : S(t+1) = S(t)*[1 + mu*dt + chol(gamma)*dWt]
      DO i = 1,N
c
c        step on each factor
         CALL GBMSTEP ( p, mudt, chol, sqrdt, dwork(pdw), dwork(pdstep))
c
c        computing price(n+1) = price(n)*(step+1)
         ks = pdstep
         DO k = 1,p
            trails(i + 1,k) = trails(i,k) * ( dwork(ks) + 1 )
            ks = ks + 1
         ENDDO
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine  MCFSIM                                    
c
c     One step of Monte-Carlo simulation by factorial method
c
c------------------------------------------------------------------------
      SUBROUTINE MCFSIM ( USERPO, P, Q, alpha, beta, vol,
     &                    painit, pfinit, T, N,
     &                    sqrdt, mudt, chol,
     &                    iusdat, dusdat, dwork, iuwork, duwork,
     &                    price, vpof, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       USERPO    : external payoff (entry pointer)            procedure
c       P         : number of asset(s) (P > 0)                   integer
c       Q         : number of factor(s) (Q > 0)                  integer
c       alpha[P]  : alpha coefficient(s) (P)                      double
c       beta[Q,P] : beta coefficient(s) (Q*P)                     double
c       vol[P]    : regression error(s) volatility (P)            double
c       painit[P] : asset(s) initial prices (P)                   double
c       pfinit[Q] : factor(s) initial prices (Q)                  double
c       H         : horizon (number of dates)                    integer
c       T         : maturity (in years)                           double
c       N         : number of steps to maturity (N > 0)          integer
c       sqrdt     : sqrt of deltat (discret. step in years)       double
c       mudt[P]   : mu(i)*dt, (P)                                 double
c       chol[P,P] : choleski factor of covariance (P*P)           double
c       iusdat    : integer user data for USERPO                 integer
c       dusdat    : double user data for USERPO                   double
c
c     WORKSPACE 
c       dwork     : (2*N + 4)*Q + (2*N + 1)*P                     double
c       iuwork    : integer USERPO workspace                     integer
c       duwork    : double USERPO workspace                       double
c
c     OUTPUT 
c       price     : payoff price                                  double
c       vpof[Q]   : factors returns (Q)                            double
c       info      : diagnostic argument                          integer
c
c     CALL   
c       MCSIM     : one Monte-Carlo simulation
c       HLRETM    : Log-returns
c       RHAPTV    : historical assets values as a function
c                  of values of factors with the APT regression
c                  and volatility
c       HLPRIM    : Log-prices
c       USERPO    : payoff function (provided by user)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function
      EXTERNAL USERPO
c
c     arguments i/o
      INTEGER P, Q, N, info
      INTEGER iusdat(*)
      DOUBLE PRECISION T, sqrdt, price
      DOUBLE PRECISION alpha(*), beta(*), vol(*), vpof(*), painit(*), 
     &                 pfinit(*), mudt(*), chol(*), dusdat(*)
c
c     workspaces
      INTEGER iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER Hsim, ndsim, nrsim, i, ipfin, pdfpri, pdfret, 
     &        pdapri, pdaret, pdw
      DOUBLE PRECISION pfi, eps
      PARAMETER ( eps = 1.0E-12 )
c     INTEGER k,j      
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdfpri= 1
c     pdfpri: pointer for factor price trails, ((N+1)*Q)
      pdfret= pdfpri + ( (N + 1)*Q )
c     pdfret: pointer for factor return trails, (N*Q)
      pdaret= pdfret + ( N*Q )
c     pdaret: pointer for asset return trails, (N*P)
      pdapri= pdaret + ( N*P )
c     pdapri: pointer for asset price trails,  ((N+1)*P)
      pdw   = pdapri + ( (N + 1)*P )
c     pdw   : pointer for MCSIM workspaces ( 3*Q )
c
c     Total size of dwork array = ((N+1)*Q) + (N*Q)
c                               + (N*P) + ((N+1)*P)
c                               + (3*Q)
c                               = ((2*N+1)*Q) + ((2*N+1)*P)
c                               + (3*Q)
c                               = (2*N + 4)*Q + (2*N+1)*P
c
c-----------------------------------------------------------------------
c
c     trails of the factors : dS(t) = S(t)*[mu*dt + chol(gamma)*dWt]
      CALL MCSIM ( Q, N, sqrdt, pfinit, mudt, chol, 
     &             dwork(pdw), dwork(pdfpri) )
c
c                   PROBLEM IS HERE !!!!!!
c****************************************************      
c      open(unit=1,file='simul_factor.txt',status='unknown')
c      do i=1,N*Q
c		write(1,*)  dwork(pdfpri + i - 1)
c	end do   
c	close(unit=1)
c****************************************************	     
c
c     Log-returns of factors
      Hsim = 1
      CALL HLRETM ( N+1, Q, dwork(pdfpri), Hsim, dwork(pdfret), info )
      nrsim = N + 1 - Hsim 
c
c****************************************************      
c      open(unit=1,file='factor1.txt',status='unknown')
c      do i=1,N*Q
c		write(1,*)  dwork(pdfret + i - 1)
c	end do   
c	close(unit=1)
c****************************************************	
c
c     Log-returns of assets from regression coefficients
      CALL RHAPTV ( N, P, Q, dwork(pdfret), alpha, beta,
     &              vol, dwork(pdaret) )
c
c     prices of assets from Log-returns
      CALL HLPRIM ( N, P,  dwork(pdaret), painit,
     &              ndsim, dwork(pdapri), info )
c
c     payoff evaluation
      CALL USERPO ( P, ndsim, dwork(pdapri), T,
     &              iusdat, dusdat, iuwork, duwork, price, info )
c
c     computing factors payoff
      ipfin = pdfpri + N
      DO i = 1,Q
         pfi = pfinit(i)
         IF (abs(pfi) .LT. eps) THEN
            info = -402
            RETURN
         ENDIF
         vpof(i) = ( dwork(ipfin) - pfi ) / pfi
         ipfin = ipfin + N + 1
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine  MCFSIMM                                    
c
c     Computes one simulation of Monte-Carlo with factors, multipayoff
c
c-----------------------------------------------------------------------
      SUBROUTINE MCFSIMM ( USERPOM, npotot, indpo,
     &                     nasset, nfact, alpha, beta, vol,
     &                     painit, pfinit, matur, nstep,
     &                     sqdelt, mudelt, chol,
     &                     liudat, iusdat, ldudat, dusdat,
     &                     dwork, iuwork, duwork,
     &                     price, vpof, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            USERPOM: entry point of an external subroutine
c                     multipayoff provided by the user         procedure
c            npotot : number of possible payoffs                 integer
c            indpo  : indicator of used payoffs, vectot(npotot)  integer
c                     used if indpo(i)=1
c            nasset : number of assets                           integer
c            nfact  : number of factors                          integer
c            alpha  : alphas (regression of assets by factors),
c                     vector(nasset)                              double
c            beta   : betas  (regression of assets by factors),
c                     matrix(nfact*nasset)                        double
c            vol    : volatility of errors of regression
c                     vector(nasset)                              double
c            painit : initial prices of assets, vector(nasset)    double
c            pfinit : initial prices of factors, vector(nfact)    double
c            horiz  : horizon return in number of dates          integer
c            matur  : maturity (in years)                         double
c            nstep  : number of steps to maturity (in years)     integer
c            sqdelt : sqrt of deltat (discret. step in years)     double
c            mudelt : mu(i)*deltat, vector(nelem)                 double
c            chol   : choleski factor of covariance,
c                     matrix(nasset*nasset)                       double
c            liudat : size of integer user data for USERPOM      integer
c            iusdat : user data for USERPOM,  matrix(liudat,*)   integer
c                     payoff function provided by the user
c            ldudat : size of double prec. user data for USERPOM integer
c            dusdat : user data for USERPOM, matrix(ldudat,*)     double
c                     payoff function provided by the user
c
c     WORKSPACE 
c            dwork  : vector(
c                     (2*nstep+4)*nfact) + ((2*nstep+1)*nasset) ) double
c            iuwork : user workspace for payoff function, vector integer
c            duwork : user workspace for payoff function, vector  double
c
c     OUTPUT 
c            price  : price                                       double
c            vpof   : factors payoff, vector(nfact)               double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        MCSIM   : computes one simulation of Monte-Carlo
c        HLRETM  : Computes the Horizon Log-RETurns on a Matrix
c        RHAPTV  : computes historical assets values as a function
c                  of values of factors with the APT regression
c                  and volatility
c        HLPRIM  : Computes the Horizon Log-PRIces on a Matrix
c        USERPOM : Computes the multipayoff provided by the user
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external function
      EXTERNAL USERPOM
c
c     i/o parameters
      INTEGER npotot, nasset, nfact, nstep, liudat, ldudat, info
      INTEGER indpo(*), iusdat(*)
      DOUBLE PRECISION matur, sqdelt, price
      DOUBLE PRECISION alpha(*), beta(*), vol(*), vpof(*), painit(*), 
     &                 pfinit(*), mudelt(*), chol(*), dusdat(*)
c
c     workspaces
      INTEGER iuwork(*)
      DOUBLE PRECISION dwork(*), duwork(*)
c
c     local variables
      INTEGER i, horsim, ndsim, nrsim, ipfin, pdfpri, pdfret, pdapri,
     &        pdaret, pdw
      DOUBLE PRECISION pfi, eps
      PARAMETER ( eps = 1.0E-12 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdfpri= 1
c     pdfpri: pointer for factor price trails,  matrix((nstep+1)*nfact)
      pdfret= pdfpri + ((nstep+1)*nfact)
c     pdfret: pointer for factor return trails,  matrix(nstep*nfact)
      pdaret= pdfret + (nstep*nfact)
c     pdaret: pointer for asset return trails,  matrix(nstep*nasset)
      pdapri= pdaret + (nstep*nasset)
c     pdapri: pointer for asset price trails,  matrix((nstep+1)*nasset)
      pdw   = pdapri + ((nstep+1)*nasset)
c     pdw   : pointer for local workspaces,
c             MCSIM needs (3*nfact)
c
c     Totalsize of dwork array = ((nstep+1)*nfact) + (nstep*nfact)
c                              + (nstep*nasset) + ((nstep+1)*nasset)
c                              + (3*nfact)
c                              = ((2*nstep+1)*nfact) + ((2*nstep+1)*nasset)
c                              + (3*nfact)
c
c                              = ((2*nstep+4)*nfact) + ((2*nstep+1)*nasset)
c
c-----------------------------------------------------------------------
c
c     computing trails : S(t+1) = S(t) + dSt
c                        dSt = mu*dt + chol(gamma)*dWt
      CALL MCSIM ( nfact, nstep, sqdelt, pfinit, mudelt, chol,
     &             dwork(pdw), dwork(pdfpri) )
c
c     computes Log-returns of factors
      horsim = 1
      CALL HLRETM ( nstep+1, nfact, dwork(pdfpri), horsim,
     &              dwork(pdfret), info )
      nrsim = nstep + 1 - horsim
c
c     computes Log-returns of assets from regression coefficients
      CALL RHAPTV ( nstep, nasset, nfact, dwork(pdfret), alpha, beta,
     &              vol, dwork(pdaret) )
c
c     computes prices of assets from Log-returns
      CALL HLPRIM ( nstep, nasset,  dwork(pdaret), painit,
     &              ndsim, dwork(pdapri), info )
c
c     computing assets payoff
      CALL USERPOM ( npotot, indpo,
     &               nasset, ndsim, dwork(pdapri), matur,
     &               liudat, iusdat, ldudat, dusdat, iuwork, duwork,
     &               price, info )
c
c     factors payoff
      ipfin = pdfpri + nstep
      DO i=1,nfact
         pfi = pfinit(i)
         IF (abs(pfi) .LT. eps) THEN
            info = -402
            RETURN
         ENDIF
         vpof(i) = ( dwork(ipfin) - pfi ) / pfi
         ipfin = ipfin + nstep + 1
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine MURET                                      
c
c     Diffusion parameter mu from return to maturity of the asset
c     signature to be used by DICHOT(dichotomy of a monotonous function)
c
c-----------------------------------------------------------------------
      SUBROUTINE MURET ( mu, iusdat, dusdat, ret, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c            mu     : diffusion parameter                         double
c            liudat : size of integer user data for SIMUL        integer
c            iusdat : integer user data for SIMIL                integer
c                     vector(liudat)
c                     1- nsimul   : number of simulations
c            ldudat : size of double prec. user data for SIMUL   integer
c            dusdat : double precision user data for SIMEXT       double
c                     vector(ldudat)
c                     1- valini    : initial value
c                     2- sigma     : volatility
c                     3- matur     : maturity
c
c     OUTPUT 
c            ret    : return to maturity                          double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER info
      INTEGER iusdat(*)
      DOUBLE PRECISION mu, ret
      DOUBLE PRECISION dusdat(*)
c
c     local variables
      INTEGER i, nsimul
      DOUBLE PRECISION price, valini, sigma, matur, noise, rn, coeff, 
     &                 sigsqm
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      nsimul = iusdat(1)
      valini = dusdat(1)
      sigma  = dusdat(2)
      matur  = dusdat(3)
c
      coeff = (mu - 0.5*sigma*sigma) * matur
      sigsqm = sigma * sqrt(matur)
c
      price = 0.0
c
      DO i=1,nsimul
         noise = rn()
         price = price + valini * exp( coeff + sigsqm*noise )
      ENDDO
c
c     mean price at maturity
      price = price / nsimul

c     mean annualized return to maturity
      ret = ( price - valini ) / ( valini * matur )
c
      RETURN
      END
