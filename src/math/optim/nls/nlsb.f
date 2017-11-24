c=======================================================================
c
c     subroutine  NLSB                                      
c
c     Non-Linear Least-Square Solver (with lower/upper bounds)
c
c     min || y - f(x) ||**2   s.t. binf <= x <= bsup
c      x
c
c     where f(.) is a differentiable function of x
c
c------------------------------------------------------------------------
      SUBROUTINE nlsb ( simext, mxopt, nfct, xini, obs, binf, bsup,
     &                  spadim, grid, epsbfg,
     &                  liudat, iudata, ldudat, dudata,
     &                  iwork, dwork, xopt, fxopt, gxopt, info )
c-----------------------------------------------------------------------
c     INPUT :
c            SIMEXT : entry point of an external subroutine
c                     provided by the user
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            xini   : initial solution  vector(mxopt)             double
c            obs    : observations  vector(nfct)                  double
c            binf   : inferior bounds on parameters vector(mxopt) double
c            bsup   : superior bounds on parameters vector(mxopt) double
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
c     WORKSPACE 
c            iwork  : workspace                                  integer
c                     vector( 2*mxopt + 6 + liudat )
c            dwork  : double precision  workspace                 double
c                     vector( mxopt*(mxopt+13)/2
c                           + nfct*(mxopt+3+spadim) + ldudat )
c
c     OUTPUT 
c            xopt   : optimal point (mxopt)                       double
c            fxopt  : function value at optimal point             double
c            gxopt  : gradient value at optimal point (mxopt)     double 
c            info   : = 0 successful exit                        integer
c
c     CALL   
c         YV       : copy a vector in a vector
c         YVI      : copy a vector of integers in a vector of integers
c         BFGSBOXS : quasi-Newton optimizer with BFGS method
c                    (short call version)
c         SIMNLS   : subroutine used by BFGS to compute the value of the
c                   function and its gradient
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
      EXTERNAL simext, simnls, YV, YVI, bfgsboxs
c
c     arguments i/o
      INTEGER mxopt, nfct, spadim, liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION epsbfg
      DOUBLE PRECISION binf(*), bsup(*), obs(*), xini(*), grid(*), 
     &                 dudata(*), xopt(*), fxopt, gxopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER totitr, totsim, piwbf, pnfct, pspad, pliu, pldu, pinfo, 
     &        piwu, pgrad, pdwbf, pwsnls, pwsext, pobs, pgrid, pdwu
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      fxopt = 0
      CALL YV ( mxopt, xini, gxopt )
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piwbf  = 1
c     piwbf  : pointer for BFGSBOXS who needs (2*mxopt + 1)
c              so (2*mxopt + 1) more
      pnfct  = piwbf + (2*mxopt + 1)
c     pnfct  : pointer for the number of functions, ( communication 
c              with SIMEXT, user simulator )
      pspad  = pnfct + 1
c     pspad  : pointer for the space dimension ( communication 
c              with SIMEXT, user simulator )
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user,
c              ( communication with SIMEXT, user simulator )
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of SIMNLS)
      piwu   = pinfo + 1
c     piwu   : pointer for user SIMEXT workspace, so (liudat) more
c
c     Total size of iwork array = (2*mxopt + 6) + liudat
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pgrad = 1
c     pgrad : pointer for gradient at solution, so mxopt more
c
      pdwbf = pgrad + mxopt
c     pdwbf : pointer for BFGSBOXS who needs ( mxopt*(mxopt+11)/2 )
c             so ( mxopt*(mxopt+11)/2 ) more
      pwsnls = pdwbf +  mxopt*(mxopt+11)/2
c     pwsnls: pointer for workspace of SIMNLS, so (nfct*(mxopt+2)) more
      pwsext = pwsnls + nfct*(mxopt+2)
c     pwsext: pointer for workspace ( communication with SIMNLS
c             and SIMEXT, BFGSBOXS simulator ),
c             so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs  : pointer for the vector of observations ( communication 
c             with SIMNLS, BFGSBOXS simulator ), so nfct  more
c
      pgrid = pobs +  nfct
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu   : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = mxopt + mxopt*(mxopt+11)/2
c                                 + ( nfct*(mxopt+2) )
c                                 + ( nfct*(spadim+1) ) + ldudat
c          = mxopt*(mxopt+13)/2 + nfct*(mxopt+3+spadim) + ldudat
c
c-----------------------------------------------------------------------
c
c     communication with SIMNLS (BFGSBOXS simulator)
c                    and SIMEXT (user simulator)
      iwork(pnfct) = nfct
      iwork(pspad) = spadim
      iwork(pliu)  = liudat
      iwork(pldu)  = ldudat
      iwork(pinfo) = 0
      CALL YV ( nfct, obs, dwork(pobs) )
      CALL YV ( nfct*spadim, grid, dwork(pgrid) )
      CALL YV ( ldudat, dudata, dwork(pdwu) )
      CALL YVI ( liudat, iudata, iwork(piwu) )
c
c     optimization BFGS-box
      CALL bfgsboxs ( simnls, simext, mxopt, xini, epsbfg, binf, bsup,
     &                iwork(piwbf), dwork(pdwbf),
     &                xopt, fxopt, gxopt, totitr, totsim, info)
c
c     gestion of SIMNLS errors
      IF (iwork(pinfo).ne.0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END

