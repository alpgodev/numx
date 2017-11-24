c=======================================================================
c
c     subroutine RNLS                                       
c
c     Robust Non-Linear Least-Square Solver
c
c-----------------------------------------------------------------------
      SUBROUTINE rnls ( simext, mxpar, xini, nfct, moyobs, covobs,
     &                  skobs, spadim, grid, epsbfg,
     &                  liudat, iudata, ldudat, dudata,
     &                  iwork, dwork, 
     &                  xopt, fxopt, gxopt, lopt, glopt, info )
c-----------------------------------------------------------------------
c     INPUT :
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxpar  : number of parameters to find               integer
c            xini   : initial solution (mxpar)                    double
c            nfct   : number of functions (observations)         integer
c            moyobs : mean observations (nfct)                    double
c            covobs : covariance matrix of observations
c                     vectorized matrix(nfct*nfct)                double
c            skobs  : stress coefficient of observations          double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            epsbfg : BFGS precision stop test                    double
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT               integer
c                     vector(liudat)
c            ldudat : size of double prec. user data for SIMEXT  integer
c            dudata : double precision user data for SIMEXT       double
c                     vector(ldudat)
c
c     WORKSPACE :
c            iwork  : workspace                                  integer
c                     vector( 2*mxopt + 6 + 12*nfct + liudat )
c            dwork  : double precision  workspace                 double
c                     vector(  mxopt*(mxopt+19)/2
c                            + nfct*(7*nfct+mxpar+spadim+33)
c                            + 2 + ldudat )
c                     with mxopt = mxpar + 1
c
c     OUTPUT :             
c            xopt   : optimal point (mxpar)                       double
c            fxopt  : function value at optimal point             double
c            gxopt  : gradient value at optimal point (mxpar)     double 
c            lopt   : optimal Lagragian value                     double
c            glopt  : Lagrangian gradient value at optimal        double
c            info   : = 0 successful exit                        integer
c
c     CALL   :
c        IVX     : initialization at a scalar of a vector
c        INIRNLS : initial computing of RNLS
c        YV      : copy a vector in a vector
c        YVI     : copy a vector of integers in a vector of integers
c        PMV     : computing M*V = vector
c                  ( M matrix(n*m), V vector(m), gives M*V vector(n) )
c        BFGSBOXS: quasi-Newton optimizer with BFGS method
c                  ( short call version )
c        SIMRNLS : subroutine used by BFGS to compute the value of the
c                  function and its gradient
c
c-----------------------------------------------------------------------
c
c     SIMEXT SIGNATURE :
c
c     INPUT :
c            nfct   : number of functions                        integer
c            mxopt  : size of the vector solution                integer
c            xopt   : optimal solution (mxopt)                    double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            liudat : size of integers data provided by the user integer
c            iusdat : integers user data provided by the user
c                     vector(liudat)                             integer
c            ldudat : size of doubles data provided by the user  integer
c            dusdat : doubles user data provided by the user
c                     vector(ldudat)                              double
c
c     OUTPUT :
c            funct  : functions   vector(nfct)                    double
c            jacob  : jacobian    matrix(nfct*mx)                 double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      implicit none
c
      external simext, simrnls, inirnls, bfgsboxs, YV, YVI, IVX
c
      integer mxpar, nfct, spadim, liudat, ldudat, info
      double precision skobs, epsbfg
      integer iudata(*)
      double precision dudata(*)
      double precision moyobs(*), xini(*), grid(*), covobs(nfct,*)
      double precision xopt(*), fxopt, gxopt(*), lopt, glopt
c
      integer iwork(*)
      double precision dwork(*)
c
c-----------------------------------------------------------------------
c     local variables
      integer mxopt, totitr, totsim
      integer piwbf, pnfct, pspad, pliu, pldu, pinfo, piwu
      integer pxopt, pbinf, pbsup, pgrad, pdwbf, pwsnls, pwsext
      integer piwini, pdwini, pxini
      integer pgrid, psk, palpha, pbeta, pvp, prspec, prcov, pdwu
      double precision inilam, lplus, infini, eigmax
c
      parameter ( inilam = 10.0, lplus = 0.7, infini = 1.d20 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info  = 0
      fxopt = 0
      lopt  = 0
      glopt = 0
      mxopt = mxpar + 1
      call YV ( mxpar, xini, gxopt )
c
c     pointers for integer work space : iwork
c     ---------------------------------------
      piwini = 1
c     piwini : pointer for INIRNLS who needs (12*nfct), so 12*nfct more
      piwbf  = piwini + 12*nfct
c     piwbf  : pointer for BFGSBOXS who needs (2*mxopt + 1),
c              so (2*mxopt + 1) more
      pnfct  = piwbf + (2*mxopt + 1)
c     pnfct  : pointer for the number of functions, communication
c              with SIMEXT ( simulator of user )
      pspad = pnfct + 1
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
c     Total size of iwork array = ( 2*mxopt + 6 + 12*nfct + liudat )
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      prcov = 1
c     prcov : pointer for the reverse covariance matrix,
c             so (nfct*nfct) more
      pdwini = prcov + nfct*nfct
c     pdwini: pointer for INIRNLS who needs nfct*(2*nfct+27),
c             so nfct*(2*nfct+27) more
      pxini = pdwini + nfct*(2*nfct+27)
c     pxopt : pointer for vector solution of optmization, so mxopt more
      pxopt = pxini + mxopt
c     pxopt : pointer for vector solution of optmization, so mxopt more
      pbinf = pxopt + mxopt
c     pbinf : pointer for inferior bounds binf, so mxopt more
      pbsup = pbinf + mxopt
c     pbsup : pointer for superior bounds bsup, so mxopt more
      pgrad = pbsup + mxopt
c     pgrad : pointer for gradient at solution, so mxopt more
c
      pdwbf = pgrad + mxopt
c     pdwbf : pointer for BFGSBOXS who needs ( mxopt*(mxopt+11)/2 )
c             so ( mxopt*(mxopt+11)/2 ) more
c
      pwsnls = pdwbf + mxopt*(mxopt+11)/2
c     pwsnls: pointer for internal workspace of SIMRNLS,
c             so (nfct*(3*nfct+mxpar+3)) more
      pwsext = pwsnls + nfct*(3*nfct+mxpar+4)
c     pwsext: pointer for workspace ( communication with SIMRNLS,
c             and SIMEXT ),
c             so ( nfct*(spadim+nfct+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      psk = pgrid + (nfct*spadim)
c     psk   : pointer for the stress coefficient k ( communication
c             with SIMRNLS, BFGSBOXS simulator ), so 1 more
      palpha = psk + 1
c     palpha: pointer for the coefficient alpha ( communication
c             with SIMRNLS,  BFGSBOXS simulator ), so 1 more
      pbeta = palpha + 1
c     pbeta : pointer for the vector beta ( communication with SIMRNLS,
c              BFGSBOXS simulator ), so (nfct) more
      pvp = pbeta + nfct
c     pvp   : pointer for the eigenvestors matrix ( communication
c             with SIMRNLS  BFGSBOXS simulator ), so (nfct*nfct) more
      prspec = pvp + nfct*nfct
c     prspec: pointer for eigenvalues reverse vector ( communication
c             with SIMRNLS,  BFGSBOXS simulator ), so nfct more
      pdwu = prspec + nfct
c     pdwu   : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = nfct*nfct + nfct*(2*nfct+27)
c                               + 5*mxopt + mxopt*(mxopt+11)/2
c                               + nfct*(3*nfct+mxpar+4)
c                               + nfct*(spadim+nfct+2) + 2 + ldudat
c     = mxopt*(mxopt+21)/2 + nfct*(7*nfct+mxpar+spadim+33) + 2 + ldudat
c            with mxopt = mxpar + 1
c
c-----------------------------------------------------------------------
c
c     initial computing for SIMRNLS
      CALL inirnls ( nfct, moyobs, covobs, skobs,
     &               iwork(piwini), dwork(pdwini),
     &               dwork(palpha), dwork(pbeta), dwork(pvp),
     &               dwork(prspec), dwork(prcov), eigmax, info )
      IF ( info.lt.0 ) RETURN
c
c     initialization of xopt (elem. 1 to mxpar)
      CALL YV ( mxpar, xini, dwork(pxini) )
c
c     initialization of lambda > eigmax + epsilon (elem. mxopt=mxpar+1)
      dwork(pxini+mxpar) = inilam*eigmax
c
c     lower/upper bounds initialization
      CALL IVX ( mxopt, dwork(pbsup), infini )
      CALL IVX ( mxopt-1, dwork(pbinf), -infini )
      dwork(pbinf+mxopt-1) = eigmax * ( 1. + lplus )
c
c     communication with SIMEXT (user simulator)
      iwork(pnfct) = nfct
      iwork(pspad) = spadim
      iwork(pliu)  = liudat
      iwork(pldu)  = ldudat
      iwork(pinfo) = 0
      CALL YV ( nfct*spadim, grid, dwork(pgrid) )
      CALL YV ( ldudat, dudata, dwork(pdwu) )
      CALL YVI ( liudat, iudata, iwork(piwu) )
      dwork(psk) = skobs
c
c     optimization BFGS-box
      CALL bfgsboxs ( simrnls, simext, mxopt, dwork(pxini), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pdwbf), dwork(pxopt),
     &                fxopt, dwork(pgrad), totitr, totsim, info)
c
c     optimal point
      CALL YV ( mxpar, dwork(pxopt), xopt )
c
c     optimal Lagragian 
      lopt = dwork( pxopt + mxopt )      
c
c     gradient value at the minimum
      CALL YV ( mxpar, dwork(pgrad), gxopt )
      glopt = dwork( pgrad + mxopt )
c
c     gestion of SIMRNLS errors
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
