c=======================================================================
c
c     subroutine  NLSX                                       
c
c     Non-Linear Least-Square Solver - Expert version
c
c-----------------------------------------------------------------------
      SUBROUTINE nlsx ( simext, mxopt, xini, dxmin, df1, epsabs,
     &                  mode, maxitr, maxsim, binf, bsup,
     &                  nfct, obs, spadim, grid,
     &                  liudat, iudata, ldudat, dudata ,iwork, dwork, 
     &                  xopt, fxopt, gxopt, totitr, totsim, info )
c----------------------------------------------------------------------
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxopt  : number of parameters to find               integer
c            xini   : initial solution  vector(mxopt)             double
c            dxmin  : precision on the variables    vector(mxpar) double
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
c            binf   : inferior bounds on parameters vector(mxpar) double
c            bsup   : superior bounds on parameters vector(mxpar) double
c            nfct   : number of functions (observations)         integer
c            obs    : observations  vector(nfct)                  double
c            spadim : space dimension                            integer
c            grid   : grid     vector(nfct*spadim)                double
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT               integer
c                     vector(liudat)
c            ldudat : size of double prec. user data for SIMEXT  integer
c            dudata : double precision user data for SIMEXT       double
c                     vector(ldudat)
c
c     WORKSPACE 
c            iwork  : vector( 2*mxopt + 6 + liudat )             integer
c            dwork  : vector( mxopt*(mxopt+9)/2                  double
c                           + nfct*(mxopt+3+spadim) + ldudat )
c
c     OUTPUT 
c            xopt   : optimal point (mxopt)                       double
c            fxopt  : function value at the minimum               double
c            gxopt  : gradient value at the minimum,vector(mxpar) double
c            totitr : total of iterations                        integer
c            totsim : total of calls of simul                    integer
c            info   : = 0 successful exit                        integer
c
c     CALL   
c           YV      : copy a vector in a vector
c           YVI     : copy a vector of integers in a vector of integers
c           BFGSBOX : quasi-Newton optimizer with BFGS method
c           SIMNLS  : subroutine used by BFGS to compute the value of the
c                     function and its gradient
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
c            fxopt  : functions   vector(nfct)                    double
c            jacob  : jacobian    matrix(nfct*mxopt)              double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext, simnls, YV, YVI, bfgsbox
c
c     arguments i/o
      INTEGER mxopt, mode, maxitr, maxsim, nfct, spadim, totitr, totsim,
     &        liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION df1, epsabs, fxopt
      DOUBLE PRECISION obs(*), xini(*), dxmin(*), binf(*), bsup(*), 
     &                 xopt(*), gxopt(*), grid(*), dudata(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwbf, pnfct, pspad, pliu, pldu, pinfo, piwu, pdwbf, 
     &        pwsnls, pwsext, pobs, pgrid, pdwu
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
c     piwbf  : pointer for BFGSBOX who needs (2*mxopt + 1)
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
      pdwbf = 1
c     pdwbf : pointer for BFGSBOX who needs ( mxopt*(mxopt+9)/2 )
      pwsnls = pdwbf +  mxopt*(mxopt+9)/2
c     pwsnls: pointer for internal workspace of SIMNLS,
c             so (nfct*(mxopt+2)) more
      pwsext = pwsnls + nfct*(mxopt+2)
c     pwsext: pointer for workspace ( communication with SIMNLS
c             and SIMEXT ),
c             so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs  : pointer for the vector of observations ( communication 
c             with SIMNLS, BFGSBOX simulator ), so nfct  more
      pgrid = pobs +  nfct
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu  : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = mxopt*(mxopt+9)/2
c                                 + ( nfct*(mxopt+2) )
c                                 + ( nfct*(spadim+1) ) + ldudat
c          = mxopt*(mxopt+9)/2 + nfct*(mxopt+3+spadim) + ldudat
c
c-----------------------------------------------------------------------
c
c     communication with SIMNLS (BFGSBOX simulator)
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
      CALL bfgsbox ( simnls, simext, mxopt, xini,
     &               dxmin, df1, epsabs, mode, maxitr, maxsim,
     &               binf, bsup,
     &               iwork(piwbf), dwork(pdwbf),
     &               xopt, fxopt, gxopt, totitr, totsim, info )
c
c     gestion of SIMNLS errors
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END

