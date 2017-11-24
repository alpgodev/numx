c=======================================================================
c
c     subroutine RNLSBX                                     
c
c     Robust Non-Linear Least-Square Solver (lower/upper bounds) 
c     Expert version
c
c-----------------------------------------------------------------------
      SUBROUTINE rnlsbx ( simext, mxpar, xini, dxmin, df1, epsabs,
     &                    mode, maxitr, maxsim, binf, bsup,
     &                    nobs, moyobs, covobs, skobs, spadim, grid,
     &                    liudat, iudata, ldudat, dudata, iwork, dwork,
     &                    xopt, fxopt, gxopt, lopt, glopt,
     &                    totitr, totsim, info )
c-----------------------------------------------------------------------
c     INPUT 
c            simext : entry point of an external subroutine
c            mxpar  : number of parameters                       integer
c            xini   : initial solution + Lagrang.(mxpar+1)        double
c            dxmin  : precision param. + Lagrang.(mxpar+1)        double
c            df1    : estimation of the first decrease of the criter
c                     used if mode=1 computed if =<0 
c                     must be nearly abs( f(xvar)-f(xopt) ) )     double
c            epsabs : precision stop test                         double
c            mode   : type of initialization                     integer
c                     = 1 : nothing to initialize
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
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT               integer
c                     vector(liudat)
c            ldudat : size of double prec. user data for SIMEXT  integer
c            dudata : double precision user data for SIMEXT       double
c                     vector(ldudat)
c
c     WORKSPACE 
c            iwork  : vector( 2*mxopt + 6 + 12*nobs + liudat )   integer
c            dwork  : vector(  mxopt*(mxopt+13)/2                 double
c                            + nobs*(7*nobs+mxpar+spadim+33)
c                            + 2 + ldudat )
c                     with mxopt = mxpar + 1
c
c     OUTPUT 
c            xopt   : optimal point (mxpar)                       double
c            fxopt  : function value at the minimum               double
c            gxopt  : gradient value at the minimum (mxpar)       double 
c            lopt   : optimal Lagragian value                     double
c            glopt  : Lagrangian gradient value at optimal        double
c            totitr : total of iterations                        integer
c            totsim : total of calls (to simulator)              integer
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        INIRNLS : initial computing of RNLSX
c        YV      : copy a vector in a vector
c        YVI     : copy a vector of integers in a vector of integers
c        BFGSBOX : quasi-Newton optimizer with BFGS method
c                  Expert version
c        SIMRNLS : subroutine used by BFGS to compute the value of the
c                  function and its gradient
c
c-----------------------------------------------------------------------
c
c     SIMEXT SIGNATURE 
c
c     INPUT 
c            nfct   : number of functions                        integer
c            mxopt  : size of the vector solution                integer
c            xopt   : solution  vector(mxopt)                     double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            liudat : size of integers data provided by the user integer
c            iusdat : integers user data provided by the user
c                     vector(liudat)                             integer
c            ldudat : size of doubles data provided by the user  integer
c            dusdat : doubles user data provided by the user
c                     vector(ldudat)                              double
c
c     OUTPUT 
c            funct  : functions   vector(nfct)                    double
c            jacob  : jacobian    matrix(nfct*mx)                 double
c            info   : = 0 successful exit                        integer
c                    
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external functions
      EXTERNAL simext, simrnls, YV, YVI, inirnls, bfgsbox
c
c     arguments i/o
      INTEGER mxpar, mode, maxitr, maxsim, nobs, spadim, totitr, totsim,
     &        liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION df1, epsabs, skobs, funct
      DOUBLE PRECISION dudata(*), xini(*), dxmin(*), binf(*), bsup(*),
     &                 moyobs(*), grid(*), covobs(nobs,*),
     &                 xopt(*), fxopt, gxopt(*), lopt, glopt
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mxopt
      INTEGER piwbf, pnobs, pspad, pliu, pldu, pinfo, piwu,
     &        pdwbf, pwsnls, pwsext, piwini, pdwini,
     &        pgrid, psk, palpha, pbeta, pvp, prspec, prcov, pdwu,
     &        pxsol, pgrad
      DOUBLE PRECISION eigmax
c
c-----------------------------------------------------------------------
c     initializations
      info = 0
      fxopt = 0
      lopt = 0
      glopt = 0
      CALL YV ( mxpar, xini, xopt )
      CALL YV ( mxpar, xini, gxopt )
      mxopt = mxpar + 1
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piwini = 1
c     piwini : pointer for INIRNLS who needs (12*nobs), so 12*nobs more
      piwbf  = piwini + 12*nobs
c     piwbf  : pointer for BFGSBOX who needs (2*mxopt + 1),
c              so (2*mxopt + 1) more
      pnobs  = piwbf + (2*mxopt + 1)
c     pnobs  : pointer for the number of functions, communication
c              with SIMEXT ( simulator of user )
      pspad  = pnobs + 1
c     pspad  : pointer for the space dimension
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of simnls)
      piwu   = pinfo + 1
c     piwu   : pointer for user SIMEXT workspace, so (liudat) more
c
c     Total size of iwork array = ( 2*mxopt + 6 + 12*nobs + liudat )
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pxsol = 1 
c     pxsol : pointer for the optimal solution so (mxopt)
c
      pgrad = pxsol + ( mxopt )
c     pgrad : pointer for the gradient so (mxopt)
c       
      prcov = pgrad + ( mxopt ) 
c     prcov : pointer for the reverse covariance matrix,
c             so (nobs*nobs) more
      pdwini = prcov + nobs*nobs
c     pdwini: pointer for INIRNLS who needs nobs*(2*nobs+27),
c             so nobs*(2*nobs+27) more
c
      pdwbf = pdwini + nobs*(2*nobs+27)
c     pdwbf : pointer for BFGSBOX who needs ( mxopt*(mxopt+9)/2 )
c             so ( mxopt*(mxopt+9)/2 ) more
c
      pwsnls = pdwbf + mxopt*(mxopt+9)/2
c     pwsnls: pointer for internal workspace of SIMRNLS,
c             so (nobs*(3*nobs+mxpar+3)) more
      pwsext = pwsnls + nobs*(3*nobs+mxpar+4)
c     pwsext: pointer for workspace ( communication with SIMRNLS,
c             and SIMEXT ),
c             so ( nobs*(spadim+nobs+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nobs*spadim) more
      psk = pgrid + (nobs*spadim)
c     psk   : pointer for the stress coefficient k ( communication
c             with SIMRNLS, BFGSBOX simulator ), so 1 more
      palpha = psk + 1
c     palpha: pointer for the coefficient alpha ( communication
c             with SIMRNLS, BFGSBOX simulator ), so 1 more
      pbeta = palpha + 1
c     pbeta : pointer for the vector beta ( communication with SIMRNLS,
c             BFGSBOX simulator ), so (nobs) more
      pvp = pbeta + nobs
c     pvp   : pointer for the eigenvestors matrix ( communication
c             with SIMRNLS BFGSBOX simulator ), so (nobs*nobs) more
      prspec = pvp + nobs*nobs
c     prspec: pointer for eigenvalues reverse vector ( communication
c             with SIMRNLS, BFGSBOX simulator ), so nobs more
      pdwu = prspec + nobs
c     pdwu   : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = nobs*nobs + nobs*(2*nobs+27)
c                               + mxopt*(mxopt+9)/2
c                               + nobs*(3*nobs+mxpar+4)
c                               + nobs*(spadim+nobs+2) + 2 + ldudat
c     = mxopt*(mxopt+13)/2 + nobs*(7*nobs+mxpar+spadim+33) + 2 + ldudat
c            with mxopt = mxpar + 1
c
c-----------------------------------------------------------------------
c
c     initial computing for SIMRNLS
      CALL inirnls ( nobs, moyobs, covobs, skobs,
     &               iwork(piwini), dwork(pdwini),
     &               dwork(palpha), dwork(pbeta), dwork(pvp),
     &               dwork(prspec), dwork(prcov), eigmax, info )
      IF ( info .NE. 0 ) RETURN
c
c     initialization of solution
      CALL YV ( mxopt, xini, dwork(pxsol) )
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
      dwork(psk) = skobs
c
c     optimization BFGS-box
      CALL bfgsbox ( simrnls, simext, mxopt, dwork(pxsol),
     &               dxmin, df1, epsabs, mode, maxitr, maxsim,
     &               binf, bsup, iwork(piwbf), dwork(pdwbf),
     &               dwork(pxsol), funct, dwork(pgrad),
     &               totitr, totsim, info)
c
c     optimal point
      CALL YV ( mxpar, dwork(pxsol), xopt )
c
c     optimal Lagrangian
      lopt = dwork(pxsol + mxopt)
c
c     function and gradient value at the minimum
      fxopt = funct
      CALL YV ( mxpar, dwork(pgrad), gxopt )
      glopt = dwork(pgrad + mxopt)
c
c     error code
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
