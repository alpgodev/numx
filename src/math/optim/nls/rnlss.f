c=======================================================================
c
c
c     subroutine RNLSS                      
c
c     Robust Non-Linear Least Square optimization, smooth solution
c
c-----------------------------------------------------------------------
      SUBROUTINE rnlss ( simext, mxpar, xini, nfct, moyobs, covobs,
     &                   skobs, spadim, grid, epsbfg, musmo, matsmo,
     &                   liudat, iudata, ldudat, dudata,
     &                   iwork, dwork, xpar, info )
c-----------------------------------------------------------------------
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxpar  : number of parameters to find               integer
c            xini   : initial solution  vector(mxpar)             double
c            nfct   : number of functions (observations)         integer
c            moyobs : mean observations  vector(nfct)             double
c            covobs : covariance matrix of observations
c                     vectorized matrix(nfct*nfct)                double
c            skobs  : stress coefficient of observations          double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            epsbfg : BFGS precision stop test                    double
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
c                     vector( 2*mxopt + 4 + 12*nfct + liudat )
c            dwork  : double precision  workspace                 double
c                     vector(  mxopt*(mxopt+19)/2
c                            + (mxpar*(mxpar+2))
c                            + nfct*(7*nfct+mxpar+spadim+33)
c                            + 2 + ldudat )
c                     with mxopt = mxpar + 1
c
c     OUTPUT 
c            xpar   : solution  vector(mxpar)                     double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        INIRNLS : initial computing of RNLS
c        YV      : copy a vector in a vector
c        IVX     : initialization at a scalar of a vector
c        PMX     : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YVI     : copy a vector of integers in a vector of integers
c        BFGSBOXS: quasi-Newton optimizer with BFGS method
c                  ( short call version )
c        SIMRNLSS: subroutine used by BFGS to compute the value of the
c                  function and its gradient, smooth solution
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
c     OUTPUT :
c            funct  : functions   vector(nfct)                    double
c            jacob  : jacobian    matrix(nfct*mx)                 double
c            info   : = 0 successful exit                        integer
c                    
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     external functions
      EXTERNAL simext, simrnlss, YV, YVI, IVX, PMX, bfgsboxs, inirnls
c
c     arguments i/o
      INTEGER mxpar, nfct, spadim, liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION skobs, epsbfg, musmo
      DOUBLE PRECISION moyobs(*), xini(*), grid(*), matsmo(*),
     &                 covobs(nfct,*), dudata(*), xpar(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER mxopt, totitr, totsim
      INTEGER piwbf, pnfct, pspad, pliu, pldu, pinfo, piwu,
     &        pxopt, pbinf, pbsup, pgrad, pdwbf, pwsnls, pwsext,
     &        piwini, pdwini,
     &        pgrid, psk, palpha, pbeta, pvp, prspec, prcov, psmo, pdwu
      DOUBLE PRECISION inilam, lplus, infini, eigmax, funct
      PARAMETER ( inilam = 10.0, lplus = 0.0001, infini = 1.E20 )
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
c     piwini : pointer for INIRNLS who needs (12*nfct), so 12*nfct more
      piwbf = piwini + 12*nfct
c     piwbf  : pointer for BFGSBOXS who needs (2*mxopt + 1),
c             so (2*mxopt + 1) more
      pnfct = piwbf + (2*mxopt + 1)
c     pnfct  : pointer for the number of functions, communication
c              with SIMUL ( simulator of user ), so 1 more
      pspad = pnfct + 1
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
c     so size of iwork array = ( 2*mxopt + 6 + 12*nfct + liudat )
c
c
c     pointers for double precision work space  : dwork
c
      prcov = 1
c     prcov : pointer for the reverse covariance matrix,
c             so (nfct*nfct) more
      pdwini = prcov + nfct*nfct
c     pdwini: pointer for INIRNLS who needs nfct*(2*nfct+27),
c             so nfct*(2*nfct+27) more
      pxopt = pdwini + nfct*(2*nfct+27)
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
c     pwsnls: pointer for internal workspace of SIMRNLSS,
c             so (nfct*(3*nfct+mxpar+3)+2*mxpar) more
      pwsext = pwsnls + (nfct*(3*nfct+mxpar+4)+2*mxpar)
c     pwsext: pointer for workspace ( communication with SIMNLS,
c             and SIMEXT, BFGSBOXS simulator ),
c             so ( nfct*(spadim+nfct+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      psk = pgrid + (nfct*spadim)
c     psk   : pointer for the stress coefficient k ( communication
c             with SIMRNLSS, user simulator ), so 1 more
      palpha = psk + 1
c     palpha: pointer for the coefficient alpha ( communication
c             with SIMRNLSS, user simulator ), so 1 more
      pbeta = palpha + 1
c     pbeta : pointer for the vector beta ( communication with SIMRNLSS,
c             user simulator ), so (nfct) more
      pvp = pbeta + nfct
c     pvp   : pointer for the eigenvestors matrix ( communication
c             with SIMRNLSS user simulator ), so (nfct*nfct) more
      prspec = pvp + nfct*nfct
c     prspec: pointer for eigenvalues reverse vector ( communication
c             with SIMRNLSS, user simulator ), so nfct more
      psmo = prspec + nfct
c     psmo  : pointer for the smoothing matrix ( communication 
c             with SIMNLSS, BFGSBOXS simulator ), so (mxpar*mxpar) more
      pdwu = psmo + (mxpar*mxpar)
c     pdwu   : pointer for user SIMEXT workspace, so (ldudat) more
c
c        so size of dwork array = nfct*nfct + nfct*(2*nfct+27)
c                               + 4*mxopt + mxopt*(mxopt+11)/2
c                               + nfct*(3*nfct+mxpar+4)
c                               + (mxpar*mxpar) + 2*mxpar
c                               + nfct*(spadim+nfct+2) + 2 + ldudat
c        with mxopt = mxpar + 1
c                      = mxopt*(mxopt+19)/2 + (mxpar*(mxpar+2))
c                      + nfct*(7*nfct+mxpar+spadim+33) + 2 + ldudat
c
c.......................................................................
c
c     initial computing for SIMRNLSS
      CALL inirnls ( nfct, moyobs, covobs, skobs,
     &               iwork(piwini), dwork(pdwini),
     &               dwork(palpha), dwork(pbeta), dwork(pvp),
     &               dwork(prspec), dwork(prcov), eigmax, info )
      IF ( info .NE. 0 ) RETURN
c
c     initialization of xopt
      CALL YV ( mxpar, xini, dwork(pxopt) )
c
c     initialization of lambda > eigmax + epsilon
      dwork(pxopt+mxopt-1) = inilam*eigmax
c
c     initializations of bounds
      CALL IVX ( mxopt, dwork(pbsup), infini )
      CALL IVX ( mxopt-1, dwork(pbinf), -infini )
      dwork(pbinf+mxopt-1) = eigmax * ( 1 + lplus )
c
c     compute smoothing matrix musmo*matsmo
      CALL PMX ( mxpar, mxpar, matsmo, musmo, dwork(psmo) )
c
c     communication with SIMEXT (user simulator)
      iwork(pnfct) = nfct
      iwork(pspad) = spadim
      iwork(pliu) = liudat
      iwork(pldu) = ldudat
      iwork(pinfo)= 0
      CALL YV ( nfct*spadim, grid, dwork(pgrid) )
      CALL YV ( ldudat, dudata, dwork(pdwu) )
      CALL YVI ( liudat, iudata, iwork(piwu) )
      dwork(psk)= skobs
c
c     optimization BFGS
      CALL bfgsboxs ( simrnlss, simext, mxopt, dwork(pxopt), epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pdwbf), dwork(pxopt),
     &                funct, dwork(pgrad), totitr, totsim, info)
c
      CALL YV ( mxpar, dwork(pxopt), xpar )
c
c     error code
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
