c=======================================================================
c
c     Fichier d'utilitaires 
c
c     sous-programmes divers
c
c        DAYRAT  : Computes a kind of daily 'rate'
c                  from the yearly interest rate of risk-free asset
c        FREERA  : computing the risk-free rate
c        WZERO   : put zero at the values of a vector less than eps
c                  and increase the greatest element with residual
c        GESTERR : gestion of errors
c        GREPS   : computing groups in a sorted vector for a sensibility
c        POLYN   : computing a polynome (Horner algorithm)
c        INDMTI  : index of symmetric matrix in triangular inf. matrix
c
c=======================================================================
c
c     subroutine DAYRAT                      
c
c     Computes a kind of daily 'rate' from the yearly interest rate of
c     risk-free asset
c
c=======================================================================
c
c     INPUT :
c            taux   : yearly interest rate of risk-free asset     double
c            nbjour : effective number of days in a year        interger
c
c     OUTPUT :
c            taujou : daily rate                                  double
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine DAYRAT ( taux, nbjour, taujou, info )
c
      implicit none
c
      integer nbjour, info
      double precision taux, taujou
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      taujou = 0.0
c
      if ( nbjour.le.0 ) then
         info = -2
         return
      endif 
         
      taujou = taux / nbjour
c
      return
      end
c
c=======================================================================
c
c     subroutine FREERA
c
c     Computing the risk-free rate
c
c=======================================================================
c
c     INPUT :
c            nasset : size of the portfolio                      integer
c                     with risk-free asset as first component
c            wasset : portfolio weights vector(nasset)            double
c                     with risk-free to compute
c
c     OUTPUT :
c            wasset : portfolio weights vector(nasset)            double
c                     with risk-free asset as first component
c
c-----------------------------------------------------------------------
c
      subroutine FREERA ( nasset, wasset )
c
      implicit none
c
      integer nasset
      double precision wasset(*)
c
      integer i
      double precision som
c
c-----------------------------------------------------------------------
c
      som =  0.
      do i=2,nasset
         som = som + wasset(i)
      end do
      wasset(1) = 1. - som
c
      return
      end
c
c=======================================================================
c
c     subroutine WZERO
c
c     Put zero at the values of a vector less than eps
c     and increase the greatest element with residual
c
c=======================================================================
c
c     INPUT :
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                               double
c            eps    : if vect(i) < eps, vect(i) = 0.              double
c
c     OUTPUT :
c            vect   : modified vector(nvect)                      double
c
c-----------------------------------------------------------------------
c
      subroutine WZERO ( nvect, vect, eps )
c
      implicit none
c
      integer nvect
      double precision vect(*), eps
c
      integer i, imax
      double precision res, vmax
c
c-----------------------------------------------------------------------
c
      res =  0.
      vmax = vect(1)
      imax = 1
      do i=1,nvect
         if (vect(i).lt.eps) then
            res = res + vect(i)
            vect(i) = 0.
         else
            if (vect(i).gt.vmax) then
               vmax = vect(i)
               imax = i
            end if
         end if
      end do
      vect(imax) = vmax + res
c
      return
      end
c
c=======================================================================
c
c     subroutine GESTERR
c
c     Gestion of errors
c
c=======================================================================
c
c     INPUT :
c            cinfo  : information on path of errors          charcter*80
c            cpath  : addition of path                       charcter*10
c
c     OUTPUT :
c            cinfo  : information on path of errors          charcter*80
c-----------------------------------------------------------------------
c
      subroutine GESTERR ( cinfo, cpath )
c
      implicit none
c
      character*80 cinfo
      character*20 cpath
c
      integer i, imin
c
c-----------------------------------------------------------------------
c
      i = 1
      do while ( cinfo(i:i).ne.' ')
         i = i+1
      end do
      imin = min(i+20,80)
      if (i.le.80) cinfo(i:imin)=cpath
c
      return
      end
c      
c========================================================================
c
c     subroutine POLYN
c
c     Computing a polynome (Horner algorithm)
c
c=======================================================================
c
c     INPUT :
c            n      : size of the polynome                       integer
c            param  : parameters of the polynome  vector(n+1)     double
c                     in the order a(0),a(1),........,a(n-1),a(n)
c            x      : value of x                                  double
c
c     OUTPUT :
c            value  :                                             double
c                   a(n)*x**(n) + a(n-1)*x**(n-1) +...+ a(1)*x + a(0)
c
c-----------------------------------------------------------------------
c
      subroutine POLYN ( n, param, x, value )
c
      implicit none
c
      integer n
      double precision x, value
      double precision param(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      value = param(n+1)
      do i=n,1,-1
         value = value*x + param(i)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine INDMTI
c
c     index of symmetric matrix in triangular inferior matrix
c
c=======================================================================
c
c     INPUT :
c            nmat   : size of matrix                             integer
c            ks     : index in symmetric matrix                  integer
c
c     OUTPUT :
c            i      : index(row) in triangular inferior matrix   integer
c            j      : index(column) in triangular inf. matrix    integer
c
c-----------------------------------------------------------------------
c
      subroutine INDMTI ( nmat, ks, i, j )
c
      implicit none
c
      integer nmat, ks, i, j
c
      integer k
c
c-----------------------------------------------------------------------
c
      k = 1
      do i=1,nmat
         do j=1,i
            if (ks.eq.k) return
            k = k + 1
         end do
      end do
c
      return
      end

