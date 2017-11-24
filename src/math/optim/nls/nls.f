c=======================================================================
c
c     NLS (Non-Linear Least-Square Solver)
c
c     min || y - f(x) ||**2
c      x
c
c     where f(.) is a differentiable function of x
c
c-----------------------------------------------------------------------
      SUBROUTINE nls ( simext, mxopt, nfct, xini, obs, spadim, grid,
     &                 epsbfg, liudat, iudata, ldudat, dudata,
     &                 iwork, dwork, xopt, fxopt, gxopt, info )
c-----------------------------------------------------------------------
c     INPUT 
c            SIMEXT : entry point of an external subroutine (provided by user)
c            mxopt  : number of parameters to find               integer
c            nfct   : number of functions (observations)         integer
c            xini   : initial solution (mxopt)                    double
c            obs    : observations (nfct)                         double
c            spadim : space dimension                            integer
c            grid   : grid (nfct*spadim)                          double
c            epsbfg : BFGS precision stop test                    double
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT (liudat)      integer
c            ldudat : double prec. user data for SIMEXT          integer
c            dudata : double prec. user data for SIMEXT (ldudat)  double
c
c     WORKSPACE 
c            iwork  : 2*mxopt + 6 + liudat                       integer
c            dwork  :  mxopt*(mxopt+17)/2                         double
c                    + nfct*(mxopt+3+spadim) + ldudat 
c
c     OUTPUT 
c            xopt   : optimal point (mxopt)                       double
c            fxopt  : function value at the minimum               double
c            gxopt  : gradient value at the minimum (mxopt)       double 
c            info   : diagnostic argument                        integer
c
c     CALL   
c           YV      : copy vector
c           IVX     : initialization at a scalar of a vector
c           YVI     : copy a vector of integers in a vector of integers
c           BFGSBOXS: quasi-Newton optimizer with BFGS method
c                    (short call version)
c           SIMNLS  : subroutine used by BFGS to compute the value of the
c                     function and its gradient
c-----------------------------------------------------------------------
c
c     SIMEXT SIGNATURE
c
c     INPUT 
c            nfct   : number of functions                        integer
c            mxopt  : size of the vector solution                integer
c            xopt   : solution (mxopt)                            double
c            spadim : space dimension                            integer
c            grid   : grid (nfct*spadim)                          double
c            liudat : size of integers data provided by the user integer
c            iusdat : integer user data provided (liudat)        integer
c            ldudat : size of doubles data provided              integer
c            dusdat : doubles user data provided (ldudat)         double
c
c     OUTPUT 
c            funct  : functions (nfct)                            double
c            jacob  : jacobian (nfct*mx)                          double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext, simnls
c
c     arguments i/o
      INTEGER mxopt, nfct, spadim, liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION epsbfg
      DOUBLE PRECISION obs(*), xini(*), grid(*), dudata(*), xopt(*), 
     &                 fxopt, gxopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER totitr, totsim, piwbf, pnfct, pspad, pliu, pldu, pinfo, 
     &        piwu, pbinf, pbsup, pgrad, pdwbf, pwsnls, pwsext, pobs, 
     &        pgrid, pdwu
      DOUBLE PRECISION infini
      PARAMETER ( infini = 1.E20 )
c
c-----------------------------------------------------------------------
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
      pbinf = 1
c     pbinf : pointer for inferior bounds binf, so mxopt more
      pbsup = pbinf + mxopt
c     pbsup : pointer for superior bounds bsup, so mxopt more
      pgrad = pbsup + mxopt
c     pgrad : pointer for gradient at solution, so mxopt more
c
      pdwbf = pgrad + mxopt
c     pdwbf : pointer for BFGSBOXS who needs ( mxopt*(mxopt+11)/2 )
c             so ( mxopt*(mxopt+11)/2 ) more
      pwsnls = pdwbf +  mxopt*(mxopt+11)/2
c     pwsnls: pointer for internal workspace of SIMNLS,
c             so (nfct*(mxopt+2)) more
      pwsext = pwsnls + nfct*(mxopt+2)
c     pwsext: pointer for workspace ( communication with SIMNLS
c             and SIMEXT ),
c             so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs  : pointer for the vector of observations ( communication 
c             with SIMNLS, BFGSBOXS simulator ), so nfct  more
      pgrid = pobs +  nfct
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu  : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = 3*mxopt + mxopt*(mxopt+11)/2
c                                 + ( nfct*(mxopt+2) )
c                                 + ( nfct*(spadim+1) ) + ldudat
c          = mxopt*(mxopt+17)/2 + nfct*(mxopt+3+spadim) + ldudat
c
c-----------------------------------------------------------------------
c
c     lower/upper bounds initializations
      CALL IVX ( mxopt, dwork(pbinf), -infini )
      CALL IVX ( mxopt, dwork(pbsup),  infini )
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
      CALL bfgsboxs ( simnls, simext, mxopt, xini, epsbfg,
     &                dwork(pbinf), dwork(pbsup),
     &                iwork(piwbf), dwork(pdwbf), xopt,
     &                fxopt, dwork(pgrad), totitr, totsim, info)
c
c     function and gradient value at the minimum
      CALL YV ( mxopt, dwork(pgrad), gxopt )
c
c     error code
      IF (iwork(pinfo).ne.0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
