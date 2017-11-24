c=======================================================================
c
c
c     Fichier d'utilitaires                  
c
c     sous-programmes specifiques allocation avec "turn-over"
c
c        CTMVT   : computing constraints for ALLOCMVT
c                  Mean Variance ALLOCation (Markowitz)
c                  with Turn-over minimization
c
c        OPMVT   : computing optimization for ALLOCMVT
c                  Mean Variance ALLOCation (Markowitz)
c                  with Turn-over minimization
c
c        OPVART  : computing optimization for ALLOCVART
c                  optimal ALLOCation and VAlue at Risk
c                  with Turn-over minimization
c
c=======================================================================
c
c     subroutine CTMVT
c
c     Computing constraints for ALLOCMVT
c               Mean Variance ALLOCation (Markowitz)
c               with Turn-over minimization
c
c=======================================================================
c
c     INPUT :
c            nasset : size of portfolio                          integer
c            winit  : initial portfolio, vector(nasset)           double
c            rmean  : mean returns vector(nasset)                 double
c            wgcf   : weight of guarantied capital fund
c            rgcf   : return of guarantied capital fund
c            mu     : performance target                          double
c            tomax  : turn-over maxi                              double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix initial of constraints
c                     matrix (nasset*(neq+nin))                   double
c            bcst   : vector(neq+nin) initial of constraints      double
c
c     OUTPUT :
c            neqtot : number of equality constraints             integer
c            nintot : number of inequality constraints           integer
c            ccstot : matrix of constraints                       double
c                     matrix((nasset+nasset)*nctot)
c            bcstot : vector of constraints                       double
c                     vector( neq + nin + 2*nasset + 3 )
c            plin   : linear part vector(nasset+nasset)           double
c            
c     METHOD :
c                 composition of ccstot, bcstot, plin
c         - ccstot
c           | 1 | 1 | neq   |   nin   | nasset  |  nasset | 1 |  size
c           ___________________________________________________
c           |   |   |                 |         |         | - |
c           |   |   |                 |         |         | r |
c           |   |   |                 |         |         | m |
c           | 1 | 0 |      ccst       |    Id   |   -Id   | e |  nasset
c           |   |   |                 |         |         | a |  ( w(i) )
c           |   |   |                 |         |         | n |
c           |   |   |                 |         |         |   |
c           ___________________________________________________
c           |   |   |                 |         |         |   |
c           |   |   |                 |         |         |   |
c           |   |   |                 |         |         |   |
c           | 0 | 1 |       0         |   -Id   |   -Id   | 0 |  nasset
c           |   |   |                 |         |         |   |  ( u(i) )
c           |   |   |                 |         |         |   |
c           |   |   |                 |         |         |   |
c           ___________________________________________________
c
c         - bcstot
c           | 1 | 1 | neq   |   nin   | nasset  |  nasset | 1 |  size
c           ___________________________________________________
c           |1-w|tom|      bcst       | winit   |  -winit |-mu|   1
c           |   |   |                 |         |         |+rw|
c           ___________________________________________________
c
c         - plin
c           |          nasset       |        nasset         |  size
c           _________________________________________________
c           |            0          |          0            |   1
c           _________________________________________________
c
c-----------------------------------------------------------------------
c
      subroutine CTMVT ( nasset, winit, rmean, wgcf, rgcf,
     &                   mu, tomax,
     &                   neq, nin, ccst, bcst,
     &                   neqtot, nintot, ccstot, bcstot, plin )
c
      implicit none
c
      integer nasset, neq, nin, neqtot, nintot  
      double precision wgcf, rgcf, mu, tomax
      double precision winit(*), rmean(*), bcst(*), bcstot(*), plin(*)
      double precision ccst(nasset,*), ccstot(2*nasset,*)
c
      integer i, j, nopt, nctot, ncini
      integer pjsum, pjto, pjc, pjw, pjmw, pjr, piw, piu
c
c-----------------------------------------------------------------------
c
c     initialisations 
c
      ncini = neq + nin
      nopt = nasset + nasset
      neqtot = neq + 2
      nintot = nin + 2*nasset + 1
      nctot = neqtot + nintot
c
c     construction of the linear part
c
      do i=1,nopt
         plin(i) = 0.
      end do
c
c     pointers for construction of constraints matrix
c
      piw = 1
      piu = piw + nasset
c
      pjsum = 1
      pjto = pjsum + 1
      pjc = pjto + 1
      pjw = pjc + ncini
      pjmw = pjw + nasset
      pjr = pjmw + nasset
c
c     construction of the constraints matrix
c
      do i=1,nasset
c
         ccstot(piw+i-1,pjsum) = 1.0
         ccstot(piu+i-1,pjsum) = 0.0
c
         ccstot(piw+i-1,pjto) = 0.0
         ccstot(piu+i-1,pjto) = 1.0
c
         do j=1,ncini
            ccstot(piw+i-1,pjc+j-1) = ccst(i,j)
            ccstot(piu+i-1,pjc+j-1) = 0.0
         end do
c
         do j=1,nasset
            if (i.eq.j) then
               ccstot(piw+i-1,pjw+j-1) = 1.0
               ccstot(piu+i-1,pjw+j-1) = -1.0
               ccstot(piw+i-1,pjmw+j-1) = -1.0
               ccstot(piu+i-1,pjmw+j-1) = -1.0
            else
               ccstot(piw+i-1,pjw+j-1) = 0.0
               ccstot(piu+i-1,pjw+j-1) = 0.0
               ccstot(piw+i-1,pjmw+j-1) = 0.0
               ccstot(piu+i-1,pjmw+j-1) = 0.0
            end if
         end do
c
         ccstot(piw+i-1,pjr) = - rmean(i)
         ccstot(piu+i-1,pjr) = 0.0
c
      end do
c
c     construction of the constraints vector
c
      bcstot(pjsum) = 1.0 - wgcf
      bcstot(pjto) = 2.*tomax
c
      do j=1,ncini
         bcstot(pjc+j-1) = bcst(j)
      end do
      do j=1,nasset
         bcstot(pjw+j-1) = winit(j)
         bcstot(pjmw+j-1) = -winit(j)
      end do
      bcstot(pjr) = - mu + wgcf*rgcf
      bcstot(pjto) = 2.*tomax
c
      return
      end
c
c=======================================================================
c=======================================================================
c
c     subroutine OPMVT
c
c     Computing optimization for ALLOCMVT
c               Mean Variance ALLOCation (Markowitz)
c               with Turn-over minimization
c
c=======================================================================
c
c     INPUT :
c            nasset : number of risky assets                     integer
c            winit  : initial portfolio, vector(nasset+1)         double
c                     guarantied capital fund in first
c            cov    : covariance matrix(nasset*nasset)            double
c            rmean  : mean returns, vector(nasset)                double
c            rgcf   : return of guarantied capital fund
c            mu     : performance target                          double
c            tomax  : turn-over maxi (0<tomax<1)                  double
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix initial of constraints
c                     matrix (nasset*(neq+nin))                   double
c            bcst   : vector(neq+nin) initial of constraints      double
c            cinf   : assets inferior limit, vector (nasset)      double
c            csup   : assets superior limit, vector (nasset)      double
c
c     WORKSPACE :
c            iwork  : vector( 10*nasset + 2*nin + neq + 6 )      integer 
c            dwork  : matrix                                      double
c                     (  nasset*( 12*nasset + 2*neq + 2*nin + 36 )
c                        + 3*neq + 5*nin + 10  )
c
c     OUTPUT :
c            wopt   : optimal portfolio  vector(nasset+1)         double
c                     guarantied capital fund in first
c            info   : = 0 successful exit                        integer
c                     = 5 : bad value for tomax
c                     = 6 : max of iterations on alpha exceeded
c
c     CALL   :
c        TM      : computing the trace of a full matrix
c        IVX     : initialization at a scalar of a vector
c        YM      : copy a matrix in a matrix
c        YV      : copy a vector in a vector
c        OPMV    : computing optimization for ALLOCMV
c        NDVL1   : computing the L_1 norm of the difference of 2 vectors
c        IMX     : initialization at a scalar of a matrix
c        YMPMP   : copy a part of a vectorized matrix
c                  in a part of a vectorized matrix
c        CTMVT   : Computing constraints for ALLOCMVT
c        QUAPRO  : quadratic solver (vectorized version)
c
c-----------------------------------------------------------------------
c
      subroutine OPMVT ( nasset, winit, cov, rmean, rgcf, mu, tomax,
     &                   neq, nin, ccst, bcst, cinf, csup,
     &                   iwork, dwork, wopt, info )
c
      implicit none
c
      integer nasset, info, neq, nin
      double precision mu, tomax, rgcf
      double precision cov(*), winit(*), rmean(*), cinf(*), csup(*)
      double precision ccst(nasset,*), bcst(*), wopt(*)
c
      integer iwork(*)
      double precision dwork(*)
c
      integer nopt, neqtot, nintot, nctot, ncini, neqop, ncop
      integer pdbop, pdcop
      integer pdcov, pdvopt, pdinf, pdsup, pdlin, pdlagr, pdbcs, pdccs
      integer piwo, pdwo, piwq, pdwq
      double precision tomes, wgcf, mugcf
      double precision myzero, dzero, dun
c
      parameter ( myzero = 1.0d-12 )
      parameter ( dzero = 0.0d0, dun = 1.0d0 )
c
c-----------------------------------------------------------------------
c
c     initialisations 
      ncini = neq + nin
      neqop = neq + 1
      ncop = neqop + nin
      nopt = nasset + nasset
      neqtot = neq + 2
      nintot = nin + 2*nasset + 1
      nctot = neqtot + nintot
c
c     --- pointers for integer work space : iwork ---
c
      piwo = 1
c     piwo  : pointer for internal workspaces of OPMV
c             needs ( 3*nasset + 2*nin + neqop + 3 )
c             with neqop = neq + 1
c        = ( 3*nasset + 2*nin + (neq+1) + 3 )
c        = ( 3*nasset + 2*nin + neq + 4 )
      piwq = 1
c     piwq  : pointer for internal workspaces of QUAPRO
c             needs ( 3*nopt + 2*nintot + neqtot + 1 ),
c             with nopt = 2*nasset
c                neqtot = neq + 1
c                nintot = nin + 2*nasset + 2
c        = 6*nasset + 2*(nin + 2*nasset + 2) + (neq + 1) + 1 )
c        = 6*nasset + 2*nin + 4*nasset + 4 + neq + 2 )
c        = ( 10*nasset + 2*nin + neq + 6 )
c
c     together OPVM and QUAPRO need (union of space) :
c        ( 10*nasset + 2*nin + neq + 6 )
c
c     size of iwork array = ( 10*nasset + 2*nin + neq + 6 )
c
c     --- pointers for double precision work space  : dwork ---
c
c     pointers for the call of OPMV
      pdbop  = 1
c     pdbop  : pointer for the constraints vector( ncop ),
      pdcop  = pdbop + ( ncop )
c     pdcop  : pointer for the constraints matrix( nasset*ncop ),
      pdwo   = pdcop + ( nasset*ncop )
c     pdwo   : pointer for internal workspace of OPMV
c              needs(  nasset*( nasset + neqop + nin + 9 ) +
c                      2*neqop + 4*nin + 4  )
c     size of dwork array need by the call of OPMV:
c                     ncop + nasset*ncop +
c                     nasset*( nasset + neqop + nin + 9 ) +
c                     2*neqop + 4*nin + 4
c
c     with ncop = neqop + nin
c
c                   = neqop + nin + nasset*( neqop + nin ) +
c                     nasset*( nasset + neqop + nin + 9 ) +
c                     2*neqop + 4*nin + 4
c                   = nasset*( nasset + 2*neqop + 2*nin + 9 ) +
c                     3*neqop + 5*nin + 4
c
c     with neqop = neq + 1
c
c                   = nasset*( nasset + 2*neq + 2*nin + 11 ) +
c                     3*neq + 5*nin + 7
c
c     pointers for the call of QUAPRO who uses a part of the same space
c     than OPMV
      pdcov  = 1
c     pdcov  : pointer for covariance matrix( nopt*nopt ) to minimize
      pdvopt = pdcov + ( nopt*nopt )
c     pdvopt : pointer for vector( nopt ) to optimize
      pdlin  = pdvopt + ( nopt )
c     pdlin  : pointer for linear part, vector( nopt )
      pdinf  = pdlin + ( nopt )
c     pdinf  : pointer for inferior limit, vector( nopt )
      pdsup  = pdinf + ( nopt )
c     pdsup  : pointer for superior limit, vector( nopt )
      pdbcs  = pdsup + ( nopt )
c     pdbcs  : pointer for the constraints vector( nctot ),
      pdccs  = pdbcs + ( nctot )
c     pdccs  : pointer for the constraints matrix( nopt*nctot ),
      pdlagr = pdccs + ( nopt*nctot )
c     pdlagr : pointer for Lagrange multipliers vector( nopt + nctot )
      pdwq   = pdlagr + ( nopt + nctot )
c     pdwq   : pointer for QUAPRO internal workspace,
c              (nopt*nopt + 6*nopt + 2*nintot),
c     size of dwork array need by the call of QUAPRO:
c                     ( nopt*nopt ) + 4*nopt + nctot +
c                     ( nopt *( nctot ) ) + ( nopt + nctot ) +
c                     (nopt*nopt + 6*nopt + 2*nintot)
c                   = 2*( nopt*nopt ) + nopt *( nctot ) +
c                     11*nopt + 2*nctot + 2*nintot
c
c     with nctot = neqtot + nintot
c
c                   = 2*( nopt*nopt ) + nopt *( neqtot + nintot ) +
c                     11*nopt + 2*(neqtot + nintot) + 2*nintot
c                   = 2*( nopt*nopt ) + nopt *( neqtot + nintot ) +
c                     11*nopt + 2*neqtot + 4*nintot
c
c     with nopt   = 2*nasset
c          neqtot = neq + 1
c          nintot = nin + 2*nasset + 2
c
c                   = 8*( nasset*nasset ) +
c                     2*nasset*( neq + 1 + nin + 2*nasset + 2 ) +
c                     22*nasset + 2*neq + 2 + 4*(nin + 2*nasset + 2)
c                   = 8*( nasset*nasset ) +
c                     4*( nasset*nasset ) +
c                     nasset*( 2*neq + 2*nin + 6 ) +
c                     22*nasset + 2*neq + 4*nin + 8*nasset + 10
c                   = nasset*( 12*nasset + 2*neq + 2*nin + 36 ) +
c                     2*neq + 4*nin + 10
c
c     together OPVM and QUAPRO need (union of space) :
c        (  nasset*( 12*nasset + 2*neq + 2*nin + 36 ) +
c                     3*neq + 5*nin + 10  )
c
c     size of dwork array =
c        (  nasset*( 12*nasset + 2*neq + 2*nin + 36 ) +
c           3*neq + 5*nin + 10  )
c
c.......................................................................
c
c     initial computing
c
c.......................................................................
c
c     computing guarantied capital fund
c
      wgcf = winit(1)
      wopt(1) = wgcf
      mugcf = mu-wgcf*rgcf
c
c.......................................................................
c
c     optimization without turn-over constraint
c
c.......................................................................
c
c     construction of the equal constraint : SUM(w(i)) = 1-wgcf
c
      call IVX ( nasset, dwork(pdcop), dun )
      call YM ( nasset, ncini, ccst, dwork(pdcop+nasset) )
      dwork(pdbop) = 1. - wgcf
      call YV ( ncini, bcst, dwork(pdbop+1) )
c
c     optimization
      call OPMV ( nasset, cov, rmean, mugcf,
     &            neqop, nin, dwork(pdcop), dwork(pdbop), cinf, csup,
     &            iwork(piwo), dwork(pdwo), wopt(2), info )
      if(info.ne.0) return
c
c     computing turn-over
      call NDVL1 ( nasset+1, winit, wopt, tomes )
c
c     optimization with turn-over constraint
      if (2.*tomax.lt.tomes) then
c
c        construction of covariance matrix
         call IMX ( nopt, nopt, dwork(pdcov), dzero )
         call YMPMP ( nasset, nasset, nasset, nasset, cov, 1, 1,
     &                nopt, nopt, dwork(pdcov), 1, 1, info )
c        info = 0 by construction
c
c        construction of the bounds and initial vector
         call IVX ( nopt, dwork(pdinf), dzero )
         call IVX ( nopt, dwork(pdsup), dun )
         call IVX ( nopt, dwork(pdvopt), dzero )
         call YV ( nasset, cinf, dwork(pdinf) )
         call YV ( nasset, csup, dwork(pdsup) )
c
c        construction of the constraints matrix and vector,
c        and initial linear part
         call CTMVT ( nasset, winit(2), rmean, wgcf, rgcf,
     &                mu, tomax,
     &                neq, nin, ccst, bcst, neqtot, nintot,
     &                dwork(pdccs), dwork(pdbcs), dwork(pdlin) )
c
c        quadratic solver QUAPRO
         call QUAPRO ( nopt, dwork(pdcov), dwork(pdlin),
     &                 dwork(pdccs), dwork(pdbcs),
     &                 dwork(pdinf), dwork(pdsup),
     &                 neqtot, nintot, iwork(piwq), dwork(pdwq),
     &                 dwork(pdlagr), dwork(pdvopt), info )
         if(info.ne.0) return
c
c        vector solution
         call YV ( nasset, dwork(pdvopt), wopt(2) )
c
      end if
c
      return
      end
c
c=======================================================================
c=======================================================================
c
c     subroutine OPVART
c
c     Computing optimization for ALLOCVART
c               optimal ALLOCation and VAlue at Risk
c               with Turn-over minimization
c
c=======================================================================
c=======================================================================
c
c     INPUT :
c            nasset : number of assets                           integer
c            winit  : initial portfolio, vector(nasset+1)         double
c                     guarantied capital fund in first
c            cov    : covariance matrix (nasset*nasset)           double
c            rmean  : mean returns vector (nasset)                double
c            rgcf   : return of guarantied capital fund
c            neq    : number of initial equality constraints     integer
c            nin    : number of initial inequality constraints   integer
c            ccst   : matrix of constraints
c                     matrix(nasset*(neq+nin))                    double
c            bcst   : vector of constraints  vector(neq+nin)      double
c            cinf   : assets inferior limit  vector (nasset)      double
c            csup   : assets superior limit  vector (nasset)      double
c            tomax  : turn-over maxi (0<tomax<1)                  double
c            vrimax : maximum value at risk                       double
c            conflp : confidence level parameter ( 0<conflp<1 )   double
c            maxdic : maximum of dichotomy iterations            integer
c            epsdic : precision of dichotomy                      double
c
c     WORKSPACE :
c            iwork  : vector ( 10*nasset + 2*nin + neq + 6 )     integer 
c            dwork  : matrix                                      double
c                    (  nasset*(12*nasset+2*neq+2*nin+38) +
c                       3*neq + 5*nin + 11  )
c
c     OUTPUT :
c            wopt   : optimal portfolio  vector(nasset+1)         double
c                     guarantied capital fund in first
c            vriopt : optimal value at risk                       double
c            info   : = 0 successful exit                        integer
c                     = 5    : bad value for tomax
c                     = -1   : no solution for vriopt < vrimax
c                     = -5   : not ( 0<conflp<1 )
c                     = -12  : no solution with these constraints
c                     =  1   : solution non-optimal, reach maxdic
c
c     CALL   :
c        EVBOR   : computing the minimum and maximum of a vector
c        OPMVT   : computing optimization for ALLOCMVT
c        WZERO   : put zero at the values of a vector less than eps
c                  and increase the greatest element with residual
c        YV      : copy a vector in a vector
c        XV      : computing scalar product of two vectors = V1'*V2
c        OVTMCV  : computing V'*M*V = scalar
c                  ( M vectorized square matrix(n*n), V vector(n) )
c
c-----------------------------------------------------------------------
c
      subroutine OPVART ( nasset, winit, cov, rmean, rgcf,
     &                    neq, nin, ccst, bcst, cinf, csup, tomax,
     &                    vrimax, conflp, maxdic, epsdic,
     &                    iwork, dwork,
     &                    wopt, vriopt, info )
c
      implicit none
c
      integer nasset, neq, nin, maxdic, info
c
      double precision rgcf, tomax, vrimax, conflp, epsdic
      double precision vriopt
      double precision cov(*), winit(*), rmean(*), cinf(*), csup(*)
      double precision ccst(*), bcst(*), wopt(*)
c
      integer iwork(*)
      double precision dwork(*)
c
      integer niter, infmax, infqua, i
      integer piwo, pdwopt, pdmean, pdwo
c
      double precision rhopt, emret, volat, varian, proba, vrisav
      double precision rhomin, rhomax, epsvri, xvri, myzero
c
      parameter ( myzero = 1.0d-12 )
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     --- pointers for integer work space : iwork ---
c
      piwo = 1
c     piwo  : pointer for internal workspaces of OPMVT
c             needs ( 10*nasset + 2*nin + neq + 6 )
c
c             so size of iwork array = ( 10*nasset + 2*nin + neq + 6 )
c
c     --- pointers for double precision work space  : dwork ---
c
      pdmean = 1
c     pdmean : pointer for rmean used to compute VaR, vector(nasset)
      pdwopt = pdmean + ( nasset )
c     pdwopt : pointer for wopt save vector(nasset+1)
      pdwo   = pdwopt + ( nasset+1 )
c     pdwo   : pointer for internal workspace of OPMVT
c               needs(  nasset*( 12*nasset + 2*neq + 2*nin + 36 ) +
c                       3*neq + 5*nin + 10  )
c
c     so size of dwork array = (nasset) + (nasset+1) +
c                nasset*(12*nasset+2*neq+2*nin+36) + 3*neq + 5*nin + 10
c
c        = nasset*(12*nasset+2*neq+2*nin+38) + 3*neq + 5*nin + 11
c
c.......................................................................
c
c     inputs control
      if ( (conflp.le.myzero) .or. (conflp.ge.(1.-myzero)) ) then
         info = -5
         return
      end if
c
c     computing rhomin and rhomax
      call EVBOR ( nasset, rmean, rhomin, rhomax )
c
c     computing rmean =<0.to compute the VaR
      do i=1,nasset
         if (rmean(i).gt.0.) then
            dwork(pdmean+i-1) = 0.
         else
            dwork(pdmean+i-1) = rmean(i)
         end if
      end do
c
c.......................................................................
c
c     research of optimal value at risk by dichotomy
c
c.......................................................................
c
c     initializations
      vrisav = 0.
      xvri = 0.0
      rhopt = rhomax
      epsvri = 2.*epsdic
      vriopt = -1.
      infmax = -1
      infqua = 0
      niter = 0
c
      do while ( ( (niter.le.maxdic).and.(epsvri.gt.epsdic) ) .and.
     &           ( (xvri.gt.vrimax).or.(niter.ne.1) ) )
c
         niter = niter + 1
c
c        optimization ( ALLOCMVT )
         call OPMVT ( nasset, winit, cov,
     &                rmean, rgcf, rhopt, tomax,
     &                neq, nin, ccst, bcst, cinf, csup,
     &                iwork(piwo), dwork(pdwo), wopt, info )
c
c        put zero at the values of a wopt less than myzero
c        and increase the great element with residual
         call WZERO ( nasset, wopt(2), myzero )
c
c        errors gestion
         if(info.ne.0) then
            if(info.ne.-12) then
               return
            else
               infqua = - 12
            end if
         else
            infqua = 0
            call YV ( nasset+1, wopt, dwork(pdwopt) )
         end if
c
c        computing variance = w'*Cov*w
         call OVTMCV ( nasset, cov, wopt(2), varian )
c
c        computing ex-ante return = w'*rmean
         call XV ( nasset, wopt(2), dwork(pdmean), emret )
c
c        computing value at risk
         volat = sqrt(varian)
         proba = volat * sqrt( conflp/(1.-conflp) ) - emret
         xvri = min(proba,1.)
         xvri = max(xvri,0.)
c
c        saving last good result
         if(infqua.eq.0) then
            vrisav = xvri
         else
            if(niter.eq.1) then
               vrisav = xvri
               call YV ( nasset+1, wopt, dwork(pdwopt) )
               infqua = 0
               info = 0
            end if
         end if
c
c        computing stop test
         epsvri = abs(xvri-vriopt)
         vriopt = xvri
c
c        dichotomy
         if (vriopt.lt.vrimax) then
            rhomin = rhopt
            infmax = 0
         else
            rhomax = rhopt
         end if
         rhopt = ( rhomax + rhomin ) / 2.
      end do
c
c.......................................................................
c
c     gestion of output
c
c.......................................................................
c
      if (niter.gt.maxdic) info = 1
      if(infmax.ne.0) info = infmax
      if(infqua.ne.0) then
         info = infqua
         call YV ( nasset+1, dwork(pdwopt), wopt )
         vriopt = vrisav
      end if
c
      return
      end
