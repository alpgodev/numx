c=======================================================================
c
c     subroutine RNLSSX                      
c
c     Robust Non-Linear Least Square optimization, smooth solution
c
c------------------------------------------------------------------------
      SUBROUTINE rnlssx ( simext, mxpar, xini, dxmin, df1, epsabs,
     &                    mode, maxitr, maxsim, binf, bsup,
     &                    nobs, moyobs, covobs, skobs, spadim, grid,
     &                    musmo, matsmo,
     &                    liudat, iudata, ldudat, dudata, iwork, dwork,
     &                    xsol, funct, grad,
     &                    totitr, totsim, info )
c------------------------------------------------------------------------
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxpar  : number of parameters to find               integer
c            xini   : initial solution +Lagrang. vector(mxpar+1)  double
c            dxmin  : precision on  parameters +Lagrang.,
c                     vector(mxpar+1)                             double
c            df1    : estimation of the first decrease of the criter
c                     used if mode=1 computed if =<0 
c                     must be nearly abs( f(xvar)-f(xopt) ) )     double
c            epsabs : precision stop test                         double
c            mode   : type of initialization                     integer
c                     = 1 : nothing to initialize
c                     = else : see n2qn1 documentation
c            maxitr : maximum of iterations authorized           integer
c            maxsim : maximum of calls of simul authorized       integer
c                     3*maxitr if largely good
c            binf   : inferior bounds on parameters +Lagrang.,
c                     vector(mxpar+1)                             double
c            bsup   : superior bounds on parameters +Lagrang.,
c                     vector(mxpar+1)                             double
c            nobs   : number of functions (observations)         integer
c            moyobs : mean observations  vector(nobs)             double
c            covobs : covariance matrix of observations
c                     vectorized matrix(nobs*nobs)                double
c            skobs  : stress coefficient of observations          double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nobs*spadim)                double
c            musmo  : smoothing coefficient                       double 
c            matsmo : smoothing matrix  vector(mxpar*mxpar)       double 
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT               integer
c                     vector(liudat)
c            ldudat : size of double prec. user data for SIMEXT  integer
c            dudata : double precision user data for SIMEXT       double
c                     vector(ldudat)
c
c     WORKSPACE 
c            iwork  : workspace                                  integer
c                     vector( 2*mxopt + 4 + 12*nobs + liudat )
c            dwork  : double precision  workspace                 double
c                     vector(  mxopt*(mxopt+9)/2
c                            + (mxpar*(mxpar+2))
c                            + nobs*(7*nobs+mxpar+spadim+33)
c                            + 2 + ldudat )
c                     with mxopt = mxpar + 1
c
c     OUTPUT 
c            xsol   : solution +Lagrang.    vector(mxpar+1)       double
c            funct  : value of the function at the minimum        double
c            grad   : gradient solution +Lagrang. at the minimum,
c                     vector(mxpar+1)                             double
c            totitr : total of iterations                        integer
c            totsim : total of calls of simul                    integer
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        INIRNLS : initial computing of RNLS
c        YV      : copy a vector in a vector
c        IVX     : initialization at a scalar of a vector
c        PMX     : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YVI     : copy a vector of integers in a vector of integers
c        BFGSBOX : quasi-Newton optimizer with BFGS method
c                  Expert version
c        SIMRNLSS: subroutine used by BFGS to compute the value of the
c                  function and its gradient, smooth solution
c
c-----------------------------------------------------------------------
c
c     SIMEXT SIGNATURE 
c
c     INPUT 
c            nobs   : number of functions                        integer
c            mxopt  : size of the vector solution                integer
c            xopt   : solution  vector(mxopt)                     double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nobs*spadim)                double
c            liudat : size of integers data provided by the user integer
c            iusdat : integers user data provided by the user
c                     vector(liudat)                             integer
c            ldudat : size of doubles data provided by the user  integer
c            dusdat : doubles user data provided by the user
c                     vector(ldudat)                              double
c
c     OUTPUT 
c            funct  : functions   vector(nobs)                    double
c            jacob  : jacobian    matrix(nobs*mx)                 double
c            info   : = 0 successful exit                        integer
c                     
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external functions
      EXTERNAL simext, simrnlss, YV, YVI, PMX, inirnls, bfgsbox
c
c     arguments i/o
      INTEGER mxpar, mode, maxitr, maxsim, nobs, spadim, totitr, totsim,
     &        liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION df1, epsabs, skobs, musmo, funct
      DOUBLE PRECISION xini(*), dxmin(*), binf(*), bsup(*), grad(*),
     &                 moyobs(*), grid(*), covobs(nobs,*), xsol(*),
     &                 dudata(*), matsmo(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mxopt
      INTEGER piwbf, pnobs, pspad, pliu, pldu, pinfo, piwu, piwini, 
     &        pdwini,
     &        pdwbf, pwsnls, pwsext, pgrid, psk, palpha, pbeta, pvp, 
     &        prspec, prcov, psmo, pdwu
      DOUBLE PRECISION eigmax
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      mxopt = mxpar + 1
c
c     pointers for integer work space  : iwork
c
      piwini = 1
c     piwini : pointer for INIRNLS who needs (12*nobs), so 12*nobs more
      piwbf = piwini + 12*nobs
c     piwbf  : pointer for BFGSBOX who needs (2*mxopt + 1),
c             so (2*mxopt + 1) more
      pnobs = piwbf + (2*mxopt + 1)
c     pnobs  : pointer for the number of functions, communication
c              with SIMUL ( simulator of user ), so 1 more
      pspad = pnobs + 1
c     pspad  : pointer for the space dimension, so 1 more
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of simnls), so 1 more
      piwu = pinfo + 1
c     piwu  : pointer for user SIMEXT workspace, so (liudat) more
c
c     so size of iwork array = ( 2*mxopt + 6 + 12*nobs + liudat )
c
c
c     pointers for double precision work space  : dwork
c
      prcov = 1
c     prcov : pointer for the reverse covariance matrix,
c             so (nobs*nobs) more
      pdwini = prcov + nobs*nobs
c     pdwini: pointer for INIRNLS who needs nobs*(2*nobs+27)
      pdwbf =  pdwini + nobs*(2*nobs+27)
c     pdwbf : pointer for BFGSBOX who needs ( mxopt*(mxopt+9)/2 )
      pwsnls = pdwbf + mxopt*(mxopt+9)/2
c     pwsnls: pointer for internal workspace of SIMRNLSS,
c             so (nobs*(3*nobs+mxpar+3)+2*mxpar) more
      pwsext = pwsnls + (nobs*(3*nobs+mxpar+4)+2*mxpar)
c     pwsext: pointer for workspace ( communication with SIMNLS,
c             and SIMEXT, BFGSBOX simulator ),
c             so ( nobs*(spadim+nobs+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nobs*spadim) more
      psk = pgrid + (nobs*spadim)
c     psk   : pointer for the stress coefficient k ( communication
c             with SIMRNLSS, user simulator ), so 1 more
      palpha = psk + 1
c     palpha: pointer for the coefficient alpha ( communication
c             with SIMRNLSS, user simulator ), so 1 more
      pbeta = palpha + 1
c     pbeta : pointer for the vector beta ( communication with SIMRNLSS,
c             user simulator ), so (nobs) more
      pvp = pbeta + nobs
c     pvp   : pointer for the eigenvestors matrix ( communication
c             with SIMRNLSS user simulator ), so (nobs*nobs) more
      prspec = pvp + nobs*nobs
c     prspec: pointer for eigenvalues reverse vector ( communication
c             with SIMRNLSS, user simulator ), so nobs more
      psmo = prspec + nobs
c     psmo  : pointer for the smoothing matrix ( communication 
c             with SIMNLSS, BFGSBOX simulator ), so (mxpar*mxpar) more
      pdwu = psmo + (mxpar*mxpar)
c     pdwu   : pointer for user SIMEXT workspace, so (ldudat) more
c
c        so size of dwork array = nobs*nobs + nobs*(2*nobs+27)
c                               + mxopt*(mxopt+9)/2
c                               + nobs*(3*nobs+mxpar+4)
c                               + (mxpar*mxpar) + 2*mxpar
c                               + nobs*(spadim+nobs+2) + 2 + ldudat
c        with mxopt = mxpar + 1
c                      = mxopt*(mxopt+11)/2 + (mxpar*(mxpar+2))
c                      + nobs*(7*nobs+mxpar+spadim+33) + 2 + ldudat
c
c-----------------------------------------------------------------------
c
c     initial computing for SIMRNLSS
      CALL inirnls ( nobs, moyobs, covobs, skobs,
     &               iwork(piwini), dwork(pdwini),
     &               dwork(palpha), dwork(pbeta), dwork(pvp),
     &               dwork(prspec), dwork(prcov), eigmax, info )
      IF ( info .NE. 0 ) RETURN
c
c     initialization of solution
      CALL YV ( mxopt, xini, xsol )
c
c     compute smoothing matrix musmo*matsmo
      CALL PMX ( mxpar, mxpar, matsmo, musmo, dwork(psmo) )
c
c     communication with SIMEXT (user simulator)
      iwork(pnobs) = nobs
      iwork(pspad) = spadim
      iwork(pliu) = liudat
      iwork(pldu) = ldudat
      iwork(pinfo)= 0
      CALL YV ( nobs*spadim, grid, dwork(pgrid) )
      CALL YV ( ldudat, dudata, dwork(pdwu) )
      CALL YVI ( liudat, iudata, iwork(piwu) )
      dwork(psk)= skobs
c
c     optimization BFGS
      CALL bfgsbox ( simrnlss, simext, mxopt, xsol,
     &               dxmin, df1, epsabs, mode, maxitr, maxsim,
     &               binf, bsup, iwork(piwbf), dwork(pdwbf),
     &               xsol, funct, grad, totitr, totsim, info)
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
