c=======================================================================
c
c     Utilities Functions
c
c     check SDLS constraints integrity
c
c-----------------------------------------------------------------------
c
c        VCSDLS   : verification de l'integrite des contraintes
c                   de la methode SDLS
c
c        VCSDLSCE : verification de l'integrite des contraintes
c                   de la methode SDLSCE
c
c        VCSDLSCI : verification de l'integrite des contraintes
c                   de la methode SDLSCI
c
c        VCSDLSCOR: verification de l'integrite des contraintes
c                   de la methode SDLSCOR
c
c        VCSDLSG  : verification de l'integrite des contraintes
c                   de la methode SDLSG
c
c=======================================================================
c
c     subroutine  VCSDLS
c
c     verification de l'integrite des contraintes de la methode SDLS
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            nmat   : dimension of the matrix                    integer
c            mct    : number of constraints                      integer
c            aict   : mct Ai symmetric matrices(nmat*nmat)
c                     of constraints  vector(mct*nmat*nmat)       double
c            bct    : vector(mct) of constraints                  double
c            xmat   : symmetric matrix solution
c            eps    : precision test                              double
c
c     OUTPUT :
c            nvcte  : number of non valid equal constraints      integer
c            indnve : non valid equal constr., vector(nvcte)     integer
c            info   : = 0 successful exit                        integer
c                     = -1 bad integrity of constraints
c
c     CALL   :
c        XM      : computing scalar product of two matrices
c                  (M1 matrix n*m, M2 matrix n*m, gives M1.M2= scalar)
c
c-----------------------------------------------------------------------
c
      subroutine vcsdls ( nmat, mct, aict, bct, xmat, eps,
     &                    nvcte, indnve, info )
c
      implicit none
      integer nmat, mct, nvcte, info
      integer indnve(*)
      double precision eps
      double precision aict(*), bct(*), xmat(nmat,*)
      integer sm, k, ka
      double precision scapro, difeq
c
      sm = nmat*nmat
      info = 0
      nvcte = 0
      ka = 1
      do k=1,mct
         call XM ( nmat, nmat, aict(ka), xmat, scapro )
         difeq = abs( scapro - bct(k) )
         if (difeq.gt.eps) then
            info = -1
            nvcte = nvcte + 1
            indnve(nvcte) = k
         end if
         ka = ka + sm
      end do
      return
      end
c
c=======================================================================
c
c     subroutine  VCSDLSCE
c
c     verification de l'integrite des contraintes de la methode SDLSCE
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            imatce : symmetric matrix(nmat*nmat)                integer
c                     equal constraints indicators
c                     if imatce(i,j) = 1 : xmat(i,j) = cmat(i,j)
c                     if imatce(i,j) = 2 : xmat(i,j) = matce(i,j)
c                     only the lower part is used
c            matce  : symmetric matrix(nmat*nmat)                 double
c                     equal constraint value if imatce(i,j) = 2
c                     only the lower part is used
c            xmat   : symmetric matrix solution
c            eps    : precision test                              double
c
c     OUTPUT :
c            nvcte  : number of non valid equal constraints      integer
c            indnve : indices of non valid equal constr.,
c                     vector(2,nvcte)                            integer
c            info   : = 0 successful exit                        integer
c                     = -1 bad integrity of constraints
c
c-----------------------------------------------------------------------
c
      subroutine vcsdlsce ( nmat, cmat, imatce, matce, xmat, eps,
     &                      nvcte, indnve, info )
c
      implicit none
      integer nmat, nvcte, info
      integer imatce(nmat,*), indnve(2,*)
      double precision eps
      double precision cmat(nmat,*), matce(nmat,*), xmat(nmat,*)
      integer i, j, ice
      double precision difeq
c
c     initialization
      info = 0
      nvcte = 0
c
      do i=1,nmat
         do j=1,i
            ice = imatce(i,j)
c
            if (ice.eq.1) then
               difeq = abs( cmat(i,j) - xmat(i,j) )
               if (difeq.gt.eps) then
                  info = -1
                  nvcte = nvcte + 1
                  indnve(1,nvcte) = i
                  indnve(2,nvcte) = j
               end if
            end if
c
            if (ice.eq.2) then
               difeq = abs( matce(i,j) - xmat(i,j) )
               if (difeq.gt.eps) then
                  info = -1
                  nvcte = nvcte + 1
                  indnve(1,nvcte) = i
                  indnve(2,nvcte) = j
               end if
            end if
c
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine  VCSDLSCI
c
c     verification de l'integrite des contraintes de la methode SDLSCI
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            nmat   : dimension of the matrix                    integer
c            cmat   : symmetric matrix(nmat*nmat) to optimize     double
c            imatc  : symmetric matrix(nmat*nmat)                integer
c                     constraints indicators (see method)
c                     only the lower part is used
c            matc   : symmetric matrix(nmat*nmat)                 double
c                     equal constraint value if imatc(i,j) = 2
c                     inequal constraint tolerance if imatc(i,j) = 3
c                     inequal constraint up value if imatc(i,j) = 4
c                     only the lower part is used
c            matlci : symmetric matrix(nmat*nmat)                 double
c                     inequal constraint low value if imatc(i,j) = 4
c                     only the lower part is used
c            xmat   : symmetric matrix(nmat*nmat) solution        double
c            eps    : precision test                              double
c
c     OUTPUT :
c            nvcte  : number of non valid equal constraints      integer
c            indnve : non valid equal constr., matrix(2,nvcte)   integer
c            nvcti  : number of non valid inequal constraints    integer
c            indnvi : non valid inequal constr., matrix(2,nvcti) integer
c            info   : = 0 successful exit                        integer
c                     = -1 bad integrity of constraints
c
c     METHOD :
c            if imatc(i,j) = 1, xmat(i,j) = cmat(i,j)
c            if imatc(i,j) = 2, xmat(i,j) = matc(i,j)
c            if imatc(i,j) = 3,
c               cmat(i,j)-matc(i,j) < xmat(i,j) < cmat(i,j)+matc(i,j)
c            if imatc(i,j) = 4, matlci(i,j) < xmat(i,j) < matc(i,j)
c
c-----------------------------------------------------------------------
c
      subroutine vcsdlsci ( nmat, cmat, imatc, matc, matlci, xmat, eps,
     &                      nvcte, indnve, nvcti, indnvi, info )
c
      implicit none
      integer nmat, nvcte, nvcti, info
      integer imatc(nmat,*), indnve(2,*), indnvi(2,*)
      double precision eps
      double precision cmat(nmat,*), matc(nmat,*), matlci(nmat,*)
      double precision xmat(nmat,*)
      integer i, j, ice
      double precision difeq, diflow, difup
c
c     initialization
      info = 0
      nvcte = 0
      nvcti = 0
c
      do i=1,nmat
         do j=1,i
            ice = imatc(i,j)
c
            if (ice.eq.1) then
               difeq = abs( cmat(i,j) - xmat(i,j) )
               if (difeq.gt.eps) then
                  info = -1
                  nvcte = nvcte + 1
                  indnve(1,nvcte) = i
                  indnve(2,nvcte) = j
               end if
            end if
c
            if (ice.eq.2) then
               difeq = abs( matc(i,j) - xmat(i,j) )
               if (difeq.gt.eps) then
                  info = -1
                  nvcte = nvcte + 1
                  indnve(1,nvcte) = i
                  indnve(2,nvcte) = j
               end if
            end if
c
            if (ice.eq.3) then
               diflow = (cmat(i,j)-matc(i,j)) - xmat(i,j)
               if (diflow.gt.eps) then
                  info = -1
                  nvcti = nvcti + 1
                  indnvi(1,nvcti) = -i
                  indnvi(2,nvcti) = j
               end if
               difup = xmat(i,j) - (cmat(i,j)+matc(i,j))
               if (difup.gt.eps) then
                  info = -1
                  nvcti = nvcti + 1
                  indnvi(1,nvcti) = i
                  indnvi(2,nvcti) = j
               end if
            end if
c
            if (ice.eq.4) then
               diflow = matlci(i,j) - xmat(i,j)
               if (diflow.gt.eps) then
                  info = -1
                  nvcti = nvcti + 1
                  indnvi(1,nvcti) = -i
                  indnvi(2,nvcti) = j
               end if
               difup = xmat(i,j) - matc(i,j)
               if (difup.gt.eps) then
                  info = -1
                  nvcti = nvcti + 1
                  indnvi(1,nvcti) = i
                  indnvi(2,nvcti) = j
               end if
            end if
c
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine  VCSDLSCOR                                 
c
c     verification de l'integrite des contraintes de la methode SDLS
c
c-----------------------------------------------------------------------
c
c     INPUT :
c            nmat   : dimension of the matrix                    integer
c            xmat   : symmetric matrix solution
c            eps    : precision test                              double
c
c     OUTPUT :
c            nvcte  : number of non valid equal constraints      integer
c            indnve : non valid equal constr., vector(nvcte)     integer
c            info   : = 0 successful exit                        integer
c                     = -1 bad integrity of constraints
c
c-----------------------------------------------------------------------
c
      subroutine vcsdlscor ( nmat, xmat, eps, nvcte, indnve, info )
c
      implicit none
      integer nmat, nvcte, info
      integer indnve(*)
      double precision eps
      double precision xmat(nmat,*)
      integer i
      double precision difeq
c
c     initialization
      info = 0
      nvcte = 0
c
      do i=1,nmat
         difeq = abs( xmat(i,i) - 1.0 )
         if (difeq.gt.eps) then
            info = -1
            nvcte = nvcte + 1
            indnve(nvcte) = i
         end if
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine  VCSDLSG                                    
c
c     verification de l'integrite des contraintes de la methode SDLSG
c
c-----------------------------------------------------------------------
      SUBROUTINE vcsdlsg ( nmat, mcte, acte, bcte,
     &                     mcti, acti, lcti, ucti, xmat, eps,
     &                     nvcte, indnve, nvcti, indnvi, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            nmat   : dimension of the matrix                    integer
c            mcte   : number of equal constraints                integer
c            acte   : mcte Ai symmetric matrices(nmat*nmat) of
c                     equal constraints  vector(mcte*nmat*nmat)   double
c            bcte   : vector(mcte) of equal constraints           double
c            mcti   : number of inequal constraints              integer
c            acti   : mcti Bi symmetric matrices(nmat*nmat) of
c                     inequal constraints  vector(mcti*nmat*nmat) double
c            lcti   : lower bounds of inequal constraints
c                     vector(mcti)                                double
c            ucti   : upper bounds of inequal constraints
c                     vector(mcti)                                double
c            xmat   : symmetric matrix solution
c            eps    : precision test                              double
c
c     OUTPUT :
c            nvcte  : number of non valid equal constraints      integer
c            indnve : non valid equal constr., vector(nvcte)     integer
c            nvcti  : number of non valid inequal constraints    integer
c            indnvi : non valid inequal constr., vector(nvcti)   integer
c            info   : = 0 successful exit                        integer
c                     = -1 bad integrity of constraints
c
c     CALL   :
c        XM      : computing scalar product of two matrices
c                  (M1 matrix n*m, M2 matrix n*m, gives M1.M2= scalar)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER nmat, mcte, mcti, nvcte, nvcti, info
      INTEGER indnve(*), indnvi(*)
      DOUBLE PRECISION eps
      DOUBLE PRECISION acte(*), bcte(*), acti(*), lcti(*), ucti(*),
     &                 xmat(nmat,*)
c
c     local variables
      INTEGER sm, k, ka
      DOUBLE PRECISION scapro, difeq, diflow, difup
c
c     initialization
      sm = nmat*nmat
      info = 0
      nvcte = 0
      nvcti = 0
c
      ka = 1
      do k=1,mcte
         call XM ( nmat, nmat, acte(ka), xmat, scapro )
         difeq = abs( scapro - bcte(k) )
         if (difeq.gt.eps) then
            info = -1
            nvcte = nvcte + 1
            indnve(nvcte) = k
         endif
         ka = ka + sm
      enddo
c
      ka = 1
      do k=1,mcti
         call XM ( nmat, nmat, acti(ka), xmat, scapro )
         diflow = lcti(k) - scapro
         if (diflow.gt.eps) then
            info = -1
            nvcti = nvcti + 1
            indnvi(nvcti) = -k
         endif
         difup = scapro - ucti(k)
         if (difup.gt.eps) then
            info = -1
            nvcti = nvcti + 1
            indnvi(nvcti) = k
         endif
         ka = ka + sm
      enddo
c
      return
      end

