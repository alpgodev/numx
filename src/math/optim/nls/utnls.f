c=======================================================================
c
c     NLS/RNLS Utility Functions
c
c-----------------------------------------------------------------------
c
c        SIMNLS  : subroutine used by BFGS to compute the value of the
c                  function and its gradient for NLS method
c
c        SIMNLSS : subroutine used by BFGS to compute the value of the
c                  function and its gradient for NLSS method,
c                  smooth solution
c
c        SIMRNLS : subroutine used by BFGS to compute the value of the
c                  function and its gradient for RNLS method
c
c        SIMRNLSS: subroutine used by BFGS to compute the value of the
c                  function and its gradient for RNLSS method,
c                  smooth solution
c
c        INIRNLS : initial computing of RNLS
c
c-----------------------------------------------------------------------
c
c     subroutine SIMNLS
c
c     Subroutine used by BFGS to compute the value of the function
c     and its gradient for NLS method
c
c-----------------------------------------------------------------------
      SUBROUTINE simnls ( indic, simext, mxopt, xsol, funct, grad,
     &                    iwork, dwork)
c-----------------------------------------------------------------------
c
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxopt  : number of variables                        integer
c            xsol   : solution  vector(mxopt) double
c
c     OUTPUT 
c            funct  : function value                              double
c            grad   : gradient value vector(mxopt)                double
c
c     WORKSPACE 
c            iwork  : vector ( 5 + liudat )                      integer
c                     -contents in input :
c                        nfct (number of functions) in iwork(1)
c                        spadim (space dimension) in iwork(2)
c                        liudat (size of integers data provided
c                                by the user) in iwork(3)
c                        ldudat (size of doubles data provided
c                                by the user) in iwork(4)
c                        user SIMEXT workspace (liudat)
c                          at iwork(6) -> iwork(5+liudat)
c                     -contents in output :
c                        info (errors of SIMEXT and SIMNLS ) in iwork(5)
c            dwork  : vector( (nfct*(mxopt+3+spadim)) + ldudat )  double
c                     with (nfct*(mxopt+2))
c                     for internal space and the rest for communication
c                     -contents in input :
c                        vector(nfct) observations
c                           at dwork((nfct*(mxopt+2))+1)
c                        vector(nfct*spadim) grid matrix
c                           at dwork((nfct*(mxopt+3))+1)
c                        user SIMEXT workspace (ldudat)
c                           at dwork((nfct*(mxopt+3+spadim))+1)
c
c     CALL   :
c        SIMEXT  : subroutine provided by the user to compute
c                  the value of the functions vector(nfct)
c                  and the jacobian matrix(nfct*mxopt)
c                  in function of xsol vector(mxopt)
c        DV      : computing the difference of 2 vectors
c        NV2     : computing the L_2 norm square of a vector
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c        PMTV    : computing M'*V = vector
c                  (M matrix(n*m), V vector(n), gives M'*V vector(m) )
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, mxopt
      DOUBLE PRECISION funct
      DOUBLE PRECISION xsol(*), grad(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER nfct, spadim, pnfct, pspad, pliu, pldu, pinfo, piwu, 
     &        pwsext, pobs, pgrid, pfsim, pjsim, pdiff, pdwu, info
      DOUBLE PRECISION deux
      PARAMETER ( deux = 2. )
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c      
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pnfct  = 1
c     pnfct  : pointer for nfct the number of functions
      pspad  = pnfct + 1
c     pspad  : pointer for the space dimension
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of SIMEXT)
      piwu = pinfo + 1
c     piwu   : pointer for user SIMEXT data, so (liudat) more
c     pitot  = piwu + (liudat)
c              
c     Total size of iwork array (5+liudat)
c
c     initializations
      nfct         = iwork(pnfct)
      spadim       = iwork(pspad)
      iwork(pinfo) = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pfsim = 1
c     pfsim  : pointer for vector of functions value provided by
c              SIMEXT, so nfct more
      pjsim = pfsim  + nfct
c     pjsim  : pointer for the jacobian value provided by
c              SIMEXT, so nfct*mxopt more
      pdiff = pjsim + nfct*mxopt
c     pdiff  : pointer for a vector(nfct), so nfct more
      pwsext = pdiff + nfct
c     pwsext: pointer for workspace ( communication with SIMNLS,
c             and SIMEXT ),
c             so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs   : pointer for the vector obs(nfct), so nfct more
      pgrid = pobs  + nfct
c     pgrid  : pointer for grid matrix need by SIMEXT,
c              so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu   : pointer for user SIMEXT data, so (ldudat) more
c     ptot = pdwu + ldudat
c            
c     Total size of dwork array = 3*nfct + (nfct*mxopt)
c                                + (nfct*spadim) + ldudat
c            = nfct*(spadim+mxopt+3) + ldudat
c
c-----------------------------------------------------------------------
c
c     objective functions value and jacobian by SIMEXT ( provided by user )
      CALL simext ( nfct, mxopt, xsol, spadim, dwork(pgrid),
     &              iwork(pliu), iwork(piwu), iwork(pldu), dwork(pdwu),
     &              dwork(pfsim), dwork(pjsim), iwork(pinfo) )
c
c     computing vector diff(i) = funct(i) - obs(i)
      CALL DV ( nfct, dwork(pfsim), dwork(pobs), dwork(pdiff) )
c
c     computing the function value = || diff || ** 2
      CALL NV2 ( nfct, dwork(pdiff), funct )
c
c     computing the gradient value = 2. * jacob(x)'*(funct-obs)
      CALL PVX ( nfct, dwork(pdiff), deux, dwork(pdiff) )
c
c     
      CALL PMTV ( nfct, mxopt, dwork(pjsim), dwork(pdiff), grad )
c
      nfct = indic
      RETURN
      END
c
c=======================================================================
c
c     subroutine SIMNLSS
c
c     Subroutine used by BFGS to compute the value of the function
c     and its gradient for NLS method, smooth solution
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxopt     : number of variables                     integer
c            xsol   : solution  vector(mxopt) double
c
c     OUTPUT 
c            funct  : function value                              double
c            grad   : gradient value vector(mxopt)                double
c
c     WORKSPACE 
c            iwork  : vector ( 5 + liudat )                      integer
c                     -contents in input :
c                        nfct (number of functions) in iwork(1)
c                        spadim (space dimension) in iwork(2)
c                        liudat (size of integers data provided
c                                by the user) in iwork(3)
c                        ldudat (size of doubles data provided
c                                by the user) in iwork(4)
c                        user SIMEXT workspace (liudat)
c                          at iwork(6) -> iwork(5+liudat)
c                     -contents in output :
c                        info (errors of SIMEXT and SIMNLS ) in iwork(5)
c            dwork  : vector(   (nfct*(mxopt+3+spadim))
c                             + (mxopt*(mxopt+2)) + ldudat   )    double
c                     with (nfct*(mxopt+2))+2*mxopt)
c                     for internal space and the rest for communication
c                     -contents in input :
c                        vector(nfct) observations
c                           at dwork((nfct*(mxopt+2))+2*mxopt+1)
c                        vector(mxopt*mxopt) smoothing matrix
c                           at dwork((nfct*(mxopt+3))+2*mxopt+1)
c                        vector(nfct*spadim) grid matrix
c                           at dwork(  (nfct*(mxopt+3))
c                                     +(mxopt*(mxopt+2))+1  )
c                        user SIMEXT workspace (ldudat)
c                           at dwork(  (nfct*(mxopt+3+spadim))
c                                      + (mxopt*(mxopt+2)) + 1  )
c
c     CALL   
c        simext  : subroutine provided by the user to compute
c                  the value of the functions vector(nfct)
c                  and the jacobian matrix(nfct*mxopt)
c                  in function of xsol vector(mxopt)
c        DV      : computing the difference of 2 vectors
c        NV2     : computing the L_2 norm square of a vector
c        PMV     : computing M*V = vector
c                  ( M matrix(n*m), V vector(m), gives M*V vector(n) )
c        XV      : computing scalar product of two vectors = V1'*V2
c        PMTV    : computing M'*V = vector
c                  (M matrix(n*m), V vector(n), gives M'*V vector(m) )
c        SV      : computing the sum of 2 vectors
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c
c-----------------------------------------------------------------------
c
      subroutine simnlss ( indic, simext, mxopt, xsol, funct, grad,
     &                     iwork, dwork)
c
      implicit none
c
      external simext
c
      integer indic, mxopt
      double precision funct
      integer iwork(*)
      double precision xsol(*), grad(*), dwork(*)
c
      integer nfct, spadim
      integer pnfct, pspad, pliu, pldu, pinfo, piwu
      integer pfsim, pjsim, pvnf, pvmx1, pvmx2
      integer pwsext, pobs, psmo, pgrid, pdwu
      double precision cqua, csmo, deux
c
      parameter ( deux = 2. )
c
      integer info
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pnfct  = 1
c     pnfct  : pointer for nfct the number of functions
      pspad  = pnfct + 1
c     pspad  : pointer for the space dimension
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of SIMEXT)
      piwu = pinfo + 1
c     piwu   : pointer for user SIMEXT data, so (liudat) more
c     pitot  = piwu + (liudat)
c
c     Total size of iwork array (5+liudat)
c
c     initializations
      nfct = iwork(pnfct)
      spadim = iwork(pspad)
      iwork(pinfo) = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pfsim = 1
c     pfsim  : pointer for vector of functions value provided by
c              SIMEXT, so nfct more
      pjsim = pfsim  + nfct
c     pjsim  : pointer for the jacobian value provided by
c              SIMEXT, so nfct*mxopt more
      pvnf = pjsim + nfct*mxopt
c     pvnf   : pointer for a vector(nfct), so nfct more
      pvmx1 = pvnf + nfct
c     pvmx1  : pointer for a vector(mxopt), so mxopt more
      pvmx2 = pvmx1 + mxopt
c     pvmx2  : pointer for a vector(mxopt), so mxopt more
      pwsext = pvmx2 + mxopt
c     pwsext : pointer for workspace ( communication with SIMNLSS,
c              and SIMEXT ),
c              so ( nfct*(spadim+1) + ldudat ) more
      pobs = pwsext
c     pobs   : pointer for the vector obs(nfct), so nfct more
      psmo = pobs + nfct
c     psmo   : pointer for the smoothing matrix ( communication 
c              with SIMNLSS ), so (mxopt*mxopt) more
      pgrid = psmo + (mxopt*mxopt)
c     pgrid  : pointer for grid matrix need by SIMEXT,
c              so (nfct*spadim) more
      pdwu = pgrid + (nfct*spadim)
c     pdwu   : pointer for user SIMEXT data, so (ldudat) more
c     ptot = pdwu + ldudat
c
c     Total size of dwork array = 3*nfct + (nfct*mxopt) + (mxopt*mxopt)
c                                + 2*mxopt + (nfct*spadim) + ldudat
c            = nfct*(spadim+mxopt+3) + mxopt*(mxopt+2) + ldudat
c
c-----------------------------------------------------------------------
c
c     computing of functions value and jacobian by SIMEXT
c                                   ( provided by user )
      call simext ( nfct, mxopt, xsol, spadim, dwork(pgrid),
     &              iwork(pliu), iwork(piwu), iwork(pldu), dwork(pdwu),
     &              dwork(pfsim), dwork(pjsim), iwork(pinfo) )
c
c     computing vector diff(i) = funct(i) - obs(i)
      call DV ( nfct, dwork(pfsim), dwork(pobs), dwork(pvnf) )
c
c     computing || diff || ** 2
      call NV2 ( nfct, dwork(pvnf), cqua )
c
c     computing product between the smoothing matrix and xsol : S*x
      call PMV ( mxopt, mxopt, dwork(psmo), xsol, dwork(pvmx1) )
c
c     computing x'*S*x
      call XV ( mxopt, xsol, dwork(pvmx1), csmo )
c
c     computing the function value = || diff || ** 2 + x'*S*x
      funct = cqua + csmo
c
c     computing the gradient value = 2. * (jacob(x)'*(funct-obs) + S*x)
      call PMTV ( nfct, mxopt, dwork(pjsim), dwork(pvnf), dwork(pvmx2) )
      call SV ( mxopt, dwork(pvmx1), dwork(pvmx2), grad )
      call PVX ( mxopt, grad, deux, grad )
c
      nfct = indic
      return
      end
c
c=======================================================================
c
c     subroutine SIMRNLS
c
c     Subroutine used by BFGS to compute the value of the function
c     and its gradient for RNLS method
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxvar  : number of variables                        integer
c            xsol   : solution  vector(mxvar)                    double
c
c     OUTPUT :
c            funct  : function value                              double
c            grad   : gradient value  vector(mxvar)               double
c
c     WORKSPACE :
c            iwork  : vector ( 5 + liudat )                      integer
c                     -contents in input :
c                        nfct (number of functions) in iwork(1)
c                        spadim (space dimension) in iwork(2)
c                        liudat (size of integers data provided
c                                by the user) in iwork(3)
c                        ldudat (size of doubles data provided
c                                by the user) in iwork(4)
c                        user SIMEXT workspace (liudat)
c                          at iwork(6) -> iwork(5+liudat)
c                     -contents in output :
c                        info (errors of SIMEXT and SIMNLS ) in iwork(5)
c            dwork  : vector                                      double
c                     ( nfct*(4*nfct+(mxmod)+spadim+6) + 2 + ldudat )
c                     with mxmod = mxvar - 1
c                     -contents in input :
c                        vector(nfct*spadim) grid matrix
c                           at dwork( ((nfct*(3*nfct+mxmod+4)) + 1 )
c                           to dwork( ((nfct*(3*nfct+mxmod+spadim+4)) )
c                        scalar k
c                           in dwork(((nfct*(3*nfct+mxmod+spadim+4))+1)
c                        scalar alpha
c                           in dwork(((nfct*(3*nfct+mxmod+spadim+4))+2)
c                        vector(nfct) beta
c                           to dwork(((nfct*(3*nfct+mxmod+spadim+5))+2)
c                        matrix(nfct*nfct) eigenvectors matrix
c                           to dwork(((nfct*(4*nfct+mxmod+spadim+5))+2)
c                        vector(nfct) eigenvalues reverse vector
c                           to dwork(((nfct*(4*nfct+mxmod+spadim+6))+2)
c                        user SIMEXT workspace (ldudat)
c                           to dwork( ((nfct*(4*nfct+mxmod+spadim+6))
c                                     +2+ldudat)
c     CALL   :
c        simext  : subroutine provided by the user to compute
c                  the value of the functions vector(nfct)
c                  and the jacobian matrix(nfct*mxvar)
c                  in function of xsol vector(mxvar)
c        OMCDMCT : computing M*D*M'
c                  ( M square matrix n*n, D vector n of diagonal matrix,
c                   gives square matrix n*n )
c        PMV     : computing M*V = vector
c                  ( M matrix(n*m), V vector(m), gives M*V vector(n) )
c        XV      : computing scalar product of two vectors = V1'*V2
c        OVTMCV  : computing V'*M*V = scalar
c                  ( M vectorized square matrix(n*n), V vector(n) )
c        NV2     : computing the L_2 norm square of a vector
c        SV      : computing the sum of 2 vectors
c        PMTV    : computing M'*V = vector
c                  (M matrix(n*m), V vector(n), gives M'*V vector(m) )
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c
c-----------------------------------------------------------------------
c
      subroutine simrnls ( indic, simext, mxvar, xsol, funct, grad,
     &                     iwork, dwork)
c
      implicit none
c
      external simext
c
      integer indic, mxvar
      double precision funct
      integer iwork(*)
      double precision xsol(*), grad(*), dwork(*)
c
      integer nfct, spadim, mxmod, i
      integer pnfct, pspad, pliu, pldu, pinfo, piwu
      integer pwsext, pgrid, psk, palpha, pbeta, pvp, prspec, pdwu
      integer pfsim, pjsim, pdfb, painv, pdainv
      integer pvnf1, pvnf2, pdwmat
      double precision lambda, scal, rspec, dad, dbd, bad, fx2
      double precision myzero, deux
c
      parameter ( myzero = 1.d-15, deux = 2. )
c
      integer info
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pnfct  = 1
c     pnfct  : pointer for nfct the number of functions
      pspad  = pnfct + 1
c     pspad  : pointer for the space dimension
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of SIMEXT)
      piwu = pinfo + 1
c     piwu   : pointer for user SIMEXT data, so (liudat) more
c     pitot  = piwu + (liudat)
c
c     Total size of iwork array (5+liudat)
c
c     initializations
      nfct = iwork(pnfct)
      spadim = iwork(pspad)
      iwork(pinfo) = 0
      mxmod = mxvar - 1
      lambda = xsol(mxvar)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pfsim = 1
c     pfsim  : pointer for vector of functions value provided by
c              SIMEXT, so nfct more
      pjsim = pfsim + nfct
c     pjsim  : pointer for the jacobian value provided by
c              SIMEXT, so nfct*mxmod more
      pdfb = pjsim + nfct*mxmod
c     pdfb   : pointer for a vector(nfct), so nfct more
      painv = pdfb + nfct
c     painv  : pointer for matrix 1/A, so (nfct*nfct) more
      pdainv = painv + nfct*nfct
c     pdainv : pointer for matrix d(1/A)/d(l), so (nfct*nfct) more
      pdwmat = pdainv + nfct*nfct
c     pdwmat : pointer for work matrix, so (nfct*nfct) more
      pvnf1 = pdwmat + nfct*nfct
c     pvnf1  : pointer for work vector 1, so (nfct) more
      pvnf2 = pvnf1 + nfct
c     pvnf2  : pointer for work vector 2, so (nfct) more
      pwsext = pvnf2 + nfct
c     pwsext : pointer for workspace ( communication with SIMRNLS,
c              and SIMEXT ),
c              so ( nfct*(spadim+nfct+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid  : pointer for grid matrix, so (nfct*spadim) more
      psk = pgrid +  (nfct*spadim)
c     psk    : pointer for the stress coefficient k, so 1 more
      palpha = psk + 1
c     palpha : pointer for the coefficient alpha, so 1 more
      pbeta = palpha + 1
c     pbeta  : pointer for the vector beta, so (nfct) more
      pvp = pbeta + nfct
c     pvp    : pointer for the eigenvestors matrix, so (nfct*nfct) more
      prspec = pvp + nfct*nfct
c     prspec : pointer for eigenvalues reverse vector, so nfct more
      pdwu = prspec + nfct
c     pdwu   : pointer for user SIMEXT data, so (ldudat) more
c     ptot = pdwu + ldudat
c
c     Total size of dwork array = 6*nfct + nfct*spadim
c                              + 2 + (nfct*mxmod)
c                              + 4*(nfct*nfct) + ldudat
c         = nfct*(4*nfct+(mxvar-1)+spadim+6) + 2 + ldudat
c
c-----------------------------------------------------------------------
c
c     computing of functions value and jacobian by SIMEXT
c                                   ( provided by user )
      call simext ( nfct, mxmod, xsol, spadim, dwork(pgrid),
     &              iwork(pliu), iwork(piwu), iwork(pldu), dwork(pdwu),
     &              dwork(pfsim), dwork(pjsim), iwork(pinfo) )
c
      do i=1,nfct
c
c        computing vector dfb(i) = funct(i) - lambda*beta(i)
         dwork(pdfb+i-1) = dwork(pfsim+i-1) - lambda*dwork(pbeta+i-1)
c
c        computing vectors of diagonal matrix
c             vls = 1/( lambda*(1/spectrum) - Id )
c        and vsls = - (1/spectrum) * 1/( lambda*(1/spectrum) - Id )**2
         rspec = dwork(prspec+i-1)
         scal = lambda*rspec - 1.
         if ( abs(scal) .gt. myzero ) then
            dwork(pvnf1+i-1) = 1./scal
            dwork(pvnf2+i-1) = - rspec / (scal*scal)
         else
            iwork(pinfo) = -2
            return
         end if
c
      end do
c
c     computing matrix 1/A = Q * 1/( lambda*(1/spectrum) - Id ) * Q'
      call OMCDMCT ( nfct, dwork(pvp), dwork(pvnf1), dwork(pdwmat),
     &               dwork(painv) )
c
c     computing matrix B = d(1/A) / d(lambda) = Q * vsls * Q'
      call OMCDMCT ( nfct, dwork(pvp), dwork(pvnf2), dwork(pdwmat),
     &               dwork(pdainv) )
c
c     computing vector aflb = (1/A)*(f(x)-lambda*beta)
      call PMV ( nfct, nfct, dwork(painv), dwork(pdfb), dwork(pvnf1) )
c
c     computing dad = (f(x)-lambda*beta)'*(1/A)*(f(x)-lambda*beta)
      call XV ( nfct, dwork(pdfb), dwork(pvnf1), dad )
c
c     computing dbd = (f(x)-lambda*beta)'*B*(f(x)-lambda*beta)
      call OVTMCV ( nfct, dwork(pdainv), dwork(pdfb), dbd )
c
c     computing || f(x) || ** 2
      call NV2 ( nfct, dwork(pfsim), fx2 )
c
c     computing the function value value
      funct = dad + fx2 + lambda*dwork(palpha)
c
c     computing ( (1/A)*(f(x)-lambda*beta) + f(x) )
      call SV ( nfct, dwork(pvnf1), dwork(pfsim), dwork(pvnf2) )
c
c     computing the gradient value ( relative at x )
c     computing 2 * (1/A)*(f(x)-lambda*beta) + f(x))
      call PMTV ( nfct, mxmod, dwork(pjsim), dwork(pvnf2), grad )
      call PVX ( mxmod, grad, deux, grad )
c
c     computing beta'*(1/A)*(f(x)-lambda*beta)
      call XV ( nfct, dwork(pbeta), dwork(pvnf1), bad )
c
c     computing the gradient value ( relative at lambda )
      grad(mxvar) = dbd - 2.*bad + dwork(palpha)
c
      nfct = indic
      return
      end
c
c=======================================================================
c
c     subroutine SIMRNLSS
c
c     Subroutine used by BFGS to compute the value of the function
c     and its gradient for RNLS method, smooth solution
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            simext : entry point of an external subroutine
c                     provided by the user
c            mxvar  : number of variables                        integer
c            xsol   : solution  vector(mxvar)                    double
c
c     OUTPUT 
c            funct  : function value                              double
c            grad   : gradient value  vector(mxvar)               double
c
c     WORKSPACE 
c            iwork  : vector ( 5 + liudat )                      integer
c                     -contents in input :
c                        nfct (number of functions) in iwork(1)
c                        spadim (space dimension) in iwork(2)
c                        liudat (size of integers data provided
c                                by the user) in iwork(3)
c                        ldudat (size of doubles data provided
c                                by the user) in iwork(4)
c                        user SIMEXT workspace (liudat)
c                          at iwork(6) -> iwork(5+liudat)
c                     -contents in output :
c                        info (errors of SIMEXT and SIMNLS ) in iwork(5)
c            dwork  : vector                                      double
c                     ( nfct*(4*nfct+(mxmod)+spadim+6)
c                       + (mxmod*(mxmod+2)) + 2 + ldudat )
c                     with mxmod = mxvar - 1
c                     with (nfct*(3*nfct+mxmod+4)+2*mxmod)
c                     for internal space and the rest for communication
c                     -contents in input :
c                        vector(nfct*spadim) grid matrix
c                           at dwork( (nfct*(3*nfct+mxmod+4))
c                                     +2*mxmod + 1 )
c                           to dwork( (nfct*(3*nfct+mxmod+spadim+4))
c                                     +2*mxmod )
c                        scalar k
c                           in dwork( (nfct*(3*nfct+mxmod+spadim+4))
c                                     +2*mxmod+1 )
c                        scalar alpha
c                           in dwork( (nfct*(3*nfct+mxmod+spadim+4))
c                                     +2*mxmod+2 )
c                        vector(nfct) beta
c                           to dwork( (nfct*(3*nfct+mxmod+spadim+5))
c                                     +2*mxmod+2 )
c                        matrix(nfct*nfct) eigenvectors matrix
c                           to dwork( (nfct*(4*nfct+mxmod+spadim+5))
c                                     +2*mxmod+2 )
c                        vector(nfct) eigenvalues reverse vector
c                           to dwork( (nfct*(4*nfct+mxmod+spadim+6))
c                                     +2*mxmod+2 )
c                        vector(mxmod*mxmod) smoothing matrix
c                           to dwork(  ((nfct*(4*nfct+mxmod+spadim+6))
c                                      +(mxmod*(mxmod+2)+2  )
c                        user SIMEXT workspace (ldudat)
c                           to dwork( ((nfct*(4*nfct+mxmod+spadim+6))
c                                     +(mxmod*(mxmod+2)+2+ldudat)
c     CALL   
c        simext  : subroutine provided by the user to compute
c                  the value of the functions vector(nfct)
c                  and the jacobian matrix(nfct*mxvar)
c                  in function of xsol vector(mxvar)
c        OMCDMCT : computing M*D*M'
c                  ( M square matrix n*n, D vector n of diagonal matrix,
c                   gives square matrix n*n )
c        PMV     : computing M*V = vector
c                  ( M matrix(n*m), V vector(m), gives M*V vector(n) )
c        XV      : computing scalar product of two vectors = V1'*V2
c        OVTMCV  : computing V'*M*V = scalar
c                  ( M vectorized square matrix(n*n), V vector(n) )
c        NV2     : computing the L_2 norm square of a vector
c        SV      : computing the sum of 2 vectors
c        PMTV    : computing M'*V = vector
c                  (M matrix(n*m), V vector(n), gives M'*V vector(m) )
c        PVX     : computing V*X = vector
c                  ( V vector(n), X scalar, gives V*X vector(n) )
c
c-----------------------------------------------------------------------
c
      subroutine simrnlss ( indic, simext, mxvar, xsol, funct, grad,
     &                      iwork, dwork)
c
      implicit none
c
      external simext
c
      integer indic, mxvar
      double precision funct
      integer iwork(*)
      double precision xsol(*), grad(*), dwork(*)
c
      integer nfct, spadim, mxmod, i
      integer pnfct, pspad, pliu, pldu, pinfo, piwu
      integer pwsext, pgrid, psk, palpha, pbeta, pvp, prspec, psmo, pdwu
      integer pfsim, pjsim, pdfb, painv, pdainv
      integer pvnf1, pvnf2, pvmx1, pvmx2, pdwmat
      double precision lambda, scal, rspec, dad, dbd, bad, fx2, myzero
      double precision csmo, deux
c
      parameter ( myzero = 1.d-15, deux = 2. )
c
      integer info
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pnfct  = 1
c     pnfct  : pointer for nfct the number of functions
      pspad  = pnfct + 1
c     pspad  : pointer for the space dimension
      pliu   = pspad + 1
c     pliu   : pointer for size of integers data provided by the user
      pldu   = pliu + 1
c     pldu   : pointer for size of doubles data provided by the user
      pinfo  = pldu + 1
c     pinfo  : pointer for info (errors of SIMEXT)
      piwu = pinfo + 1
c     piwu   : pointer for user SIMEXT data, so (liudat) more
c     pitot  = piwu + (liudat)
c
c     Total size of iwork array (5+liudat)
c
c     initializations
      nfct = iwork(pnfct)
      spadim = iwork(pspad)
      iwork(pinfo) = 0
      mxmod = mxvar - 1
      lambda = xsol(mxvar)
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pfsim = 1
c     pfsim  : pointer for vector of functions value provided by
c              SIMEXT, so nfct more
      pjsim = pfsim + nfct
c     pjsim  : pointer for the jacobian value provided by
c              SIMEXT, so nfct*mxmod more
      pdfb = pjsim + nfct*mxmod
c     pdfb   : pointer for a vector(nfct), so nfct more
      painv = pdfb + nfct
c     painv  : pointer for matrix 1/A, so (nfct*nfct) more
      pdainv = painv + nfct*nfct
c     pdainv : pointer for matrix d(1/A)/d(l), so (nfct*nfct) more
      pdwmat = pdainv + nfct*nfct
c     pdwmat : pointer for work matrix, so (nfct*nfct) more
      pvnf1 = pdwmat + nfct*nfct
c     pvnf1  : pointer for work vector 1, so (nfct) more
      pvnf2 = pvnf1 + nfct
c     pvnf2  : pointer for work vector 2, so (nfct) more
      pvmx1 = pvnf2 + nfct
c     pvmx1  : pointer for a vector(mxmod), so mxmod more
      pvmx2 = pvmx1 + mxmod
c     pvmx2  : pointer for a vector(mxmod), so mxmod more
      pwsext = pvmx2 + mxmod
c     pwsext:  pointer for workspace ( communication with SIMRNLSS,
c              and SIMEXT ),
c              so ( nfct*(spadim+nfct+2) + 2 + ldudat ) more
      pgrid = pwsext
c     pgrid  : pointer for grid matrix, so (nfct*spadim) more
      psk = pgrid +  (nfct*spadim)
c     psk    : pointer for the stress coefficient k, so 1 more
      palpha = psk + 1
c     palpha : pointer for the coefficient alpha, so 1 more
      pbeta = palpha + 1
c     pbeta  : pointer for the vector beta, so (nfct) more
      pvp = pbeta + nfct
c     pvp    : pointer for the eigenvestors matrix, so (nfct*nfct) more
      prspec = pvp + nfct*nfct
c     prspec : pointer for eigenvalues reverse vector, so nfct more
      psmo = prspec + nfct
c     psmo   : pointer for the smoothing matrix, so (mxmod*mxmod) more
      pdwu = psmo + (mxmod*mxmod)
c     pdwu   : pointer for user SIMEXT data, so (ldudat) more
c     ptot = pdwu + ldudat
c        
c     Total size of dwork array = 6*nfct + nfct*spadim
c                              + 2 + (nfct*mxmod) + 4*(nfct*nfct)
c                              + 2*mxmod + (mxmod*mxmod) + ldudat
c        = nfct*(4*nfct+(mxmod)+spadim+6) + mxmod*(mxmod+2) + 2 + ldudat
c          with mxmod = mxvar - 1
c
c-----------------------------------------------------------------------
c
c     computing of functions value and jacobian by SIMEXT
c                                   ( provided by user )
      call simext ( nfct, mxmod, xsol, spadim, dwork(pgrid),
     &              iwork(pliu), iwork(piwu), iwork(pldu), dwork(pdwu),
     &              dwork(pfsim), dwork(pjsim), iwork(pinfo) )
c
      do i=1,nfct
c
c        computing vector dfb(i) = funct(i) - lambda*beta(i)
         dwork(pdfb+i-1) = dwork(pfsim+i-1) - lambda*dwork(pbeta+i-1)
c
c        computing vectors of diagonal matrix
c             vls = 1/( lambda*(1/spectrum) - Id )
c        and vsls = - (1/spectrum) * 1/( lambda*(1/spectrum) - Id )**2
         rspec = dwork(prspec+i-1)
         scal = lambda*rspec - 1.
         if ( abs(scal) .gt. myzero ) then
            dwork(pvnf1+i-1) = 1./scal
            dwork(pvnf2+i-1) = - rspec / (scal*scal)
         else
            iwork(pinfo) = -2
            return
         end if
c
      end do
c
c     computing matrix 1/A = Q * 1/( lambda*(1/spectrum) - Id ) * Q'
      call OMCDMCT ( nfct, dwork(pvp), dwork(pvnf1), dwork(pdwmat),
     &               dwork(painv) )
c
c     computing matrix B = d(1/A) / d(lambda) = Q * vsls * Q'
      call OMCDMCT ( nfct, dwork(pvp), dwork(pvnf2), dwork(pdwmat),
     &               dwork(pdainv) )
c
c     computing vector aflb = (1/A)*(f(x)-lambda*beta)
      call PMV ( nfct, nfct, dwork(painv), dwork(pdfb), dwork(pvnf1) )
c
c     computing dad = (f(x)-lambda*beta)'*(1/A)*(f(x)-lambda*beta)
      call XV ( nfct, dwork(pdfb), dwork(pvnf1), dad )
c
c     computing dbd = (f(x)-lambda*beta)'*B*(f(x)-lambda*beta)
      call OVTMCV ( nfct, dwork(pdainv), dwork(pdfb), dbd )
c
c     computing || f(x) || ** 2
      call NV2 ( nfct, dwork(pfsim), fx2 )
c
c     computing product between the smoothing matrix and xsol : S*x
      call PMV ( mxmod, mxmod, dwork(psmo), xsol, dwork(pvmx1) )
c
c     computing x'*S*x
      call XV ( mxmod, xsol, dwork(pvmx1), csmo )
c
c     computing the function value value
      funct = dad + fx2 + lambda*dwork(palpha) + csmo
c
c     computing ( (1/A)*(f(x)-lambda*beta) + f(x) )
      call SV ( nfct, dwork(pvnf1), dwork(pfsim), dwork(pvnf2) )
c
c     computing the non smoothing half gradient value ( relative at x )
c     computing   jacob(x) * ( (1/A)*(f(x)-lambda*beta) + f(x) )
      call PMTV ( nfct, mxmod, dwork(pjsim), dwork(pvnf2),
     &            dwork(pvmx2) )
c
c     computing the gradient value ( relative at x )
c     computing 2.* ( grad + S*x )
      call SV ( mxmod, dwork(pvmx1), dwork(pvmx2), grad )
      call PVX ( mxmod, grad, deux, grad )
c
c     computing beta'*(1/A)*(f(x)-lambda*beta)
      call XV ( nfct, dwork(pbeta), dwork(pvnf1), bad )
c
c     computing the gradient value ( relative at lambda )
      grad(mxvar) = dbd - 2.*bad + dwork(palpha)
c
      nfct = indic
      return
      end
c
c=======================================================================
c
c     subroutine INIRNLS
c
c     initial computing of RNLS
c
c-----------------------------------------------------------------------
c
c     INPUT 
c            nfct   : number of functions                        integer
c            moyobs : mean observations  vector(nfct)             double
c            covobs : covariance matrix of observations
c                     vectorized matrix(nfct*nfct)                double
c            skobs  : stress coefficient of observations          double
c
c     WORKSPACE 
c            iwork  : vector (12*nfct)                           integer
c            dwork  : vector  nfct*(2*nfct+27)                    double
c
c     OUTPUT 
c            alpha  : coeficient alpha                            double
c            beta   : vector beta  vector(nfct)                   double
c            vpcov  : eigenvectors matrix
c                     vectorized matrix(nfct*nfct)                double
c            rspec  : reverse spectrum   vector(nfct)             double
c            rcov   : reverse covariance matrix
c                     vectorized matrix(nfct*nfct)                double
c            eigmax : eigenvalue max                              double
c            info   : = 0 successful exit                        integer
c
c     CALL   
c        SCHURS  : SCHUR decomposition of a symmetric matrix
c        OMCDMCT : M*D*M'
c                  ( M square matrix n*n, D vector n of diagonal matrix,
c                   gives square matrix n*n )
c        PMV     : M*V = vector
c                  ( M matrix(n*m), V vector(m), gives M*V vector(n) )
c        XV      : computing scalar product of two vectors = V1'*V2
c
c-----------------------------------------------------------------------
c
      subroutine inirnls ( nfct, moyobs, covobs, skobs, iwork, dwork,
     &                     alpha, beta, vpcov, rspec, rcov, eigmax,
     &                     info )
c
      implicit none
c
      integer nfct, info
      double precision skobs, alpha, eigmax
      integer iwork(*)
      double precision covobs(nfct,*), moyobs(*), dwork(*)
      double precision vpcov(nfct,*), rcov(nfct,*), rspec(*), beta(*)
c
      integer i
      integer pmat, peig, pdw, piw
      double precision scal, myzero
c
      parameter ( myzero = 1.d-10 )
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw  = 1
c     piw    : pointer for SCHURS who needs (12*nfct), so 12*nfct  more
c     pitot  = piw + 12*nfct
c
c     Total size of iwork array (12*nfct)
c
c     pointers for double precision work space  : work
c     ------------------------------------------------
      pmat = 1
c     pmat  : pointer for a work matrix, so (nfct*nfct) more
      peig = pmat + nfct*nfct
c     peig  : pointer for a work vector, so nfct  more
      pdw = peig + nfct
c     pdw   : pointer for SCHURS who needs nfct*(nfct+26),
c             so nfct*(nfct+26) more
c     ptot = pdw + nfct*(nfct+26)
c          
c     Total size of dwork array = nfct*(2*nfct+27)
c
c-----------------------------------------------------------------------
c
c     Schur factorization
      call SCHURS ( nfct, covobs, iwork(piw), dwork(pdw),
     &              dwork(peig), vpcov, info )
      if ( info.ne.0 ) return
c
c     computing rspec = 1/spectrum   and   lmax > max(eigenvalues)
      eigmax = 0.
      do i=1,nfct
         scal  = dwork(peig+i-1)
         eigmax = max(eigmax,scal)
         if ( abs(scal) .le. myzero ) then
            info = -2
            return
         end if
         rspec(i) = 1. / scal
      end do
c
c     computing 1/covobs
      call OMCDMCT ( nfct, vpcov, rspec, dwork(pmat), rcov )
c
c     computing vector(nfct) beta = (1/covobs)*moyobs
      call PMV ( nfct, nfct, rcov, moyobs, beta )
c
c     computing moyobs'*(1/covobs)*moyobs = moyobs'*beta
      call XV ( nfct, moyobs, beta, scal )
c
c     computing alpha = ( k**2 - moyobs'*(1/covobs)*moyobs )
      alpha = skobs*skobs - scal
c
      return
      end
