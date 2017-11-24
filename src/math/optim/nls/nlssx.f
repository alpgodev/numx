c=======================================================================
c
c     subroutine  NLSSX                                      
c
c     Non-Linear Least-Square, smooth solution - Expert version
c
c-----------------------------------------------------------------------
      SUBROUTINE nlssx ( simext, mxopt, xini, dxmin, df1, epsabs,
     &                   mode, maxitr, maxsim, binf, bsup,
     &                   nfct, obs, spadim, grid, musmo, matsmo,
     &                   liudat, iudata, ldudat, dudata ,iwork, dwork, 
     &                   xopt, funct, grad, totitr, totsim,
     &                   info )
c-----------------------------------------------------------------------
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
c            musmo  : smoothing coefficient                       double 
c            matsmo : smoothing matrix  vector(mxopt*mxopt)       double 
c            liudat : size of integer user data for SIMEXT       integer
c            iudata : integer user data for SIMEXT               integer
c                     vector(liudat)
c            ldudat : size of double prec. user data for SIMEXT  integer
c            dudata : double precision user data for SIMEXT       double
c                     vector(ldudat)
c
c     WORKSPACE 
c            iwork  : vector( 2*mxopt + 6 + liudat )             integer
c            dwork  : vector( mxopt*(3*mxopt+13)/2                double
c                           + nfct*(mxopt+3+spadim) + ldudat )
c
c     OUTPUT 
c            xopt   : solution  vector(mxopt)                     double
c            funct  : value of the function at the minimum        double
c            grad   : gradient value at the minimum,vector(mxpar) double
c            totitr : total of iterations                        integer
c            totsim : total of calls of simul                    integer
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        YV      : copy a vector in a vector
c        IVX     : initialization at a scalar of a vector
c        PMX     : computing M*X = matrix
c                  ( M matrix(n*m), X scalar, gives M*X matrix(n*m) )
c        YVI     : copy a vector of integers in a vector of integers
c        BFGSBOX : quasi-Newton optimizer with BFGS method
c                  Expert version
c        SIMNLSS : subroutine used by BFGS to compute the value of the
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
      EXTERNAL simext, simnlss, YV, PMX, YVI, bfgsbox
c
c     arguments i/o
      INTEGER mxopt, mode, maxitr, maxsim, nfct, spadim, totitr, totsim,
     &        liudat, ldudat, info
      INTEGER iudata(*)
      DOUBLE PRECISION df1, epsabs, musmo, funct
      DOUBLE PRECISION obs(*), xini(*), dxmin(*), binf(*), bsup(*),
     &                 xopt(*), grad(*), grid(*), dudata(*), matsmo(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER piwbf, pnfct, pspad, pliu, pldu, pinfo, piwu,
     &        pdwbf, pwsnls, pwsext, pobs, pgrid, psmo, pdwu
      DOUBLE PRECISION infini
      PARAMETER ( infini = 1.E20 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      CALL YV ( mxopt, xini, xopt )
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
c     pinfo  : pointer for info (errors of SIMNLSS)
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
c     pwsnls: pointer for internal workspace of SIMNLSS,
c             so (nfct*(mxopt+2)+2*mxopt) more
      pwsext = pwsnls + (nfct*(mxopt+2)+2*mxopt)
c     pwsext: pointer for workspace ( communication with SIMNLSS
c             and SIMEXT, BFGSBOX simulator ),
c             so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs  : pointer for the vector of observations ( communication 
c             with SIMNLSS, BFGSBOX simulator ), so nfct  more
      psmo = pobs + nfct
c     psmo  : pointer for the smoothing matrix ( communication 
c             with SIMNLSS, BFGSBOX simulator ), so (mxopt*mxopt) more
      pgrid = psmo + (mxopt*mxopt)
c     pgrid : pointer for grid matrix ( communication with SIMEXT,
c             user simulator ), so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu  : pointer for user SIMEXT workspace, so (ldudat) more
c
c     Total size of dwork array = mxopt*(mxopt+9)/2
c                               + ( nfct*(mxopt+2) ) + (mxopt*mxopt)
c                               + 2*mxopt
c                               + ( nfct*(spadim+1) ) + ldudat
c          = mxopt*(3*mxopt+13)/2 + nfct*(mxopt+3+spadim) + ldudat
c
c-----------------------------------------------------------------------
c
c     compute smoothing matrix musmo*matsmo
      CALL PMX ( mxopt, mxopt, matsmo, musmo, dwork(psmo) )
c
c     communication with SIMNLSS (BFGSBOX simulator)
c                    and SIMEXT (user simulator)
      iwork(pnfct) = nfct
      iwork(pspad) = spadim
      iwork(pliu)  = liudat
      iwork(pldu)  = ldudat
      iwork(pinfo) = 0
      CALL YV ( nfct, obs, dwork(pobs) )
      CALL YV ( nfct*spadim, grid, dwork(pgrid) )
c       do i=1,ldudat
c            write(*,*) 'dudata nlssx ',i,dudata(i)
c       enddo
      CALL YV ( ldudat, dudata, dwork(pdwu) )
      CALL YVI ( liudat, iudata, iwork(piwu) )
c
c     optimization BFGS
      CALL bfgsbox ( simnlss, simext, mxopt, xopt,
     &               dxmin, df1, epsabs, mode, maxitr, maxsim,
     &               binf, bsup,
     &               iwork(piwbf), dwork(pdwbf),
     &               xopt, funct, grad, totitr, totsim, info )
c
c     gestion of SIMNLSS errors
      IF (iwork(pinfo) .NE. 0) THEN
         info = iwork(pinfo)
         RETURN
      ENDIF
c
      RETURN
      END
