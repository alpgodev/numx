c=======================================================================
c
c     Fichier d'utilitaires                 
c
c     sous-programmes inpout/output pour les tests
c
c        LECRV  : lecture d'un vecteur de reels
c        LECRM  : lecture d'une matrice de reels
c        IMPRV  : imprime un vecteur de reels
c        IMPRVE : imprime un vecteur de reels (format exp)
c        IMPRM  : imprime une matrice de reels
c        IMPRME : imprime une matrice de reels (format exp)
c        IMPRMC : imprime une matrice carree de reels
c        IMPRMS : imprime une matrice vectorisee symetrique de reels
c        ARMV   : converting a matrix in vector of matrix
c        ARSMSV : converting a symmetric matrix in vector of symmetric
c                 matrix
c        ARVM   : converting a vector of matrix in matrix
c        ARSVSM : converting a vector of symmetric matrix in symmetric
c                 matrix
c        GERETU : gestion of final error in test program with return
c        GESTOP : gestion of final error in test program with stop
c        VWIN   : verification workspace - initialisation
c        VWOUT  : verification workspace - sortie
c        IMPETF : impression en-tete fonction
c        IMPETT : impression en-tete test
c
c=======================================================================
c
c     subroutine LECRV
c
c     Lecture d'un vecteur de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de lecture                           integer
c            nvect  : taille du vecteur                          integer
c
c     OUTPUT :
c            vect   : vecteur (nvect)                             double
c
c-----------------------------------------------------------------------
c
      subroutine LECRV ( iosort, nvect, vect )
c
      implicit none
c
      integer iosort, nvect
      double precision vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      read(iosort,*)(vect(i),i=1,nvect)
c
      return
      end
c
c=======================================================================
c
c     subroutine LECRM
c
c     Lecture d'une matrice de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de lecture                           integer
c            nlmat  : nombre de lignes de la matrice mat         integer
c            ncmat  : nombre de colonnes de la matrice mat       integer
c
c     OUTPUT :
c            mat    : matrice vectorisee (nlmat*ncmat)            double
c
c-----------------------------------------------------------------------
c
      subroutine LECRM ( iosort, nlmat, ncmat, mat )
c
      implicit none
c
      integer iosort, nlmat, ncmat
      double precision mat(nlmat,*)
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      do i=1,nlmat
         read(iosort,*)(mat(i,j),j=1,ncmat)
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPRV
c
c     Imprime un vecteur de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nvect  : dimension du vecteur vect
c            vect   : vecteur de dimension nvect
c
c-----------------------------------------------------------------------
c
      subroutine IMPRV ( iosort, texte, nvect, vect )
c
      implicit none
c
      integer iosort, nvect
      double precision vect(*)
      character*20 texte
c
      integer i
c
c-----------------------------------------------------------------------
c
      write(iosort,*)
      write(iosort,*) texte
      write(iosort,1) (vect(i),i=1,nvect)
      write(iosort,*)
c
    1 format(50f15.6)
c
      return
      end
c
c=======================================================================

c     subroutine IMPRVE
c
c     Imprime un vecteur de reels (format exp)
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nvect  : dimension du vecteur vect
c            vect   : vecteur de dimension nvect
c
c-----------------------------------------------------------------------
c
      subroutine IMPRVE ( iosort, texte, nvect, vect )
c
      implicit none
c
      integer iosort, nvect
      double precision vect(*)
      character*20 texte
c
      integer i
c
c-----------------------------------------------------------------------
c
      write(iosort,*)
      write(iosort,*) texte
      write(iosort,1) (vect(i),i=1,nvect)
      write(iosort,*)
c
    1 format(50e10.2)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPRM
c
c     Imprime une matrice de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nlmat  : nombre de lignes de la matrice mat
c            ncmat  : nombre de colonnes de la matrice mat
c            mat    : matrice de dimension nlmat*ncmat
c
c-----------------------------------------------------------------------
c
      subroutine IMPRM ( iosort, texte, nlmat, ncmat, mat )
c
      implicit none
c
      integer iosort, nlmat, ncmat
      double precision mat(nlmat,*)
      character*20 texte
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      write(iosort,*)
      write(iosort,*) texte
      do i=1,nlmat
         write(iosort,1)(mat(i,j),j=1,ncmat)
      end do
      write(iosort,*)
c
    1 format(50f20.15)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPRME
c
c     Imprime une matrice de reels (format exp)
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nlmat  : nombre de lignes de la matrice mat
c            ncmat  : nombre de colonnes de la matrice mat
c            mat    : matrice de dimension nlmat*ncmat
c
c-----------------------------------------------------------------------
c
      subroutine IMPRME ( iosort, texte, nlmat, ncmat, mat )
c
      implicit none
c
      integer iosort, nlmat, ncmat
      double precision mat(nlmat,*)
      character*20 texte
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      write(iosort,*)
      write(iosort,*) texte
      do i=1,nlmat
         write(iosort,1)(mat(i,j),j=1,ncmat)
      end do
      write(iosort,*)
c
    1 format(50e10.2)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPRMC
c
c     Imprime une matrice carree de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nmat   : taille de la matrice carree mat
c            mat    : matrice de dimension nmat*nmat
c
c-----------------------------------------------------------------------
c
      subroutine IMPRMC ( iosort, texte, nmat, mat )
c
      implicit none
c
      integer iosort, nmat
      double precision mat(nmat,*)
      character*20 texte
c
      integer i, j
c
c-----------------------------------------------------------------------
c
      write(iosort,*)
      write(iosort,*) texte
      do i=1,nmat
         write(iosort,1)(mat(i,j),j=1,nmat)
      end do
      write(iosort,*)
c
    1 format(50f10.4)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPRMS
c
c     Imprime une matrice symetrique vectorisee de reels
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nmat   : taille de la matrice mat
c            vmat   : matrice symetrique vectorisee (nmat*(nmat+1)/2)
c     CALL   :
c           CMSMC   : converting a vector of symmetric matrix
c                     in vector of full square matrix
c
c-----------------------------------------------------------------------
c
      subroutine IMPRMS ( iosort, texte, nmat, vmat )
c
      implicit none
c
      integer iosort, nmat
      double precision vmat(*)
      character*20 texte
c
      integer i, j
      double precision mat(nmat,nmat)
c
c-----------------------------------------------------------------------
c
      call CMSMC ( nmat, vmat, mat)
c
      write(iosort,*)
      write(iosort,*) texte
      do i=1,nmat
         write(iosort,1)(mat(i,j),j=1,nmat)
      end do
      write(iosort,*)
c
    1 format(50f10.4)
c
      return
      end
c
c=======================================================================
c
c     subroutine ARMV
c
c     Converting a matrix in vector of matrix
c
c=======================================================================
c
c     INPUT :
c            nmax   : max of first dimension of the matrix mat   integer
c            nlmat  : number of raws in matrix mat               integer
c            ncmat  : number of columns in matrix mat            integer
c            mat    : matrix (nlmat*ncmat)                        double
c
c     OUTPUT :
c            vect   : vector (nlmat*ncmat)                        double
c
c-----------------------------------------------------------------------
c
      subroutine ARMV ( nmax, nlmat, ncmat, mat, vect )
c
      implicit none
c
      integer nmax, nlmat, ncmat
      double precision mat(nmax,*), vect(*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      do j=1,ncmat
         do i=1,nlmat
            k = k + 1
            vect(k) = mat(i,j)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine ARSMSV
c
c     Converting a symmetric matrix in vector of symmetric matrix
c
c=======================================================================
c
c     INPUT :
c            nmax   : max of first dimension of the matrix mat   integer
c            nmat   : size of matrix mat                         integer
c            mat    : symetric matrix (nmat*nmat)                 double
c
c     OUTPUT :
c            nvect  : size of vector vect                        integer
c            vect   : symmetric matrix as vector(nmat*(nmat+1)/2) double
c 
c-----------------------------------------------------------------------
c
      subroutine ARSMSV ( nmax, nmat, mat, vect )
c
      implicit none
c
      integer nmax, nmat
      double precision mat(nmax,*), vect(*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      do j=1,nmat
         do i=1,j
            k = k + 1
            vect(k) = mat(i,j)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine ARVM
c
c     Converting a vector of matrix in matrix
c
c=======================================================================
c
c     INPUT :
c            vect   : vector (nlmat*ncmat)                        double
c            nmax   : max of first dimension of the matrix mat   integer
c            nlmat  : number of raws in matrix mat               integer
c            ncmat  : number of columns in matrix mat            integer
c
c     OUTPUT :
c            mat    : matrix (nlmat*ncmat)                        double
c
c-----------------------------------------------------------------------
c
      subroutine ARVM ( vect, nmax, nlmat, ncmat, mat )
c
      implicit none
c
      integer nmax, nlmat, ncmat
      double precision mat(nmax,*), vect(*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      do j=1,ncmat
         do i=1,nlmat
            k = k + 1
            mat(i,j) = vect(k)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine ARSVSM
c
c     Converting a vector of symmetric matrix in symmetric matrix
c
c=======================================================================
c
c     INPUT :
c            nmat   : size of matrix mat                         integer
c            vect   : symmetric matrix as vector(nmat*(nmat+1)/2) double
c            nmax   : max of first dimension of the matrix mat   integer
c
c     OUTPUT :
c            mat    : symmetric matrix (nmat*nmat)                double
c
c-----------------------------------------------------------------------
c
      subroutine ARSVSM ( nmat, vect, nmax, mat )
c
      implicit none
c
      integer nmax, nmat
      double precision mat(nmax,*), vect(*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      do j=1,nmat
         do i=1,j
            k = k + 1
            mat(i,j) = vect(k)
            if ( i.ne.j ) mat(j,i) = vect(k)
         end do
      end do
c
      return
      end
c
c=======================================================================
c
c     subroutine GERETU
c
c     gestion of final error in test program with return
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de sortie                            integer
c            finmes : final message                         character*20
c            cinfo  : path of error                         character*80
c            info   : error number                               integer
c
c     CALL   :
c        GESTERR : gestion of errors
c
c-----------------------------------------------------------------------
c
      subroutine GERETU ( iosort, finmes, cinfo, info )
c
      implicit none
c
      integer iosort, info
      character*20 finmes
      character*80 cinfo
c
c-----------------------------------------------------------------------
c
      if (info.ne.0) then
         call GESTERR ( cinfo, finmes )
         write(iosort,*)'**************warning********************'
         write(iosort,*)'             error ',info,' in : ',cinfo
         write(iosort,*)'**************warning********************'
         write(iosort,*)
      end if
c
      return
      end
c
c=======================================================================
c
c     subroutine GESTOP
c
c     gestion of final error in test program with stop
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de sortie                            integer
c            finmes : final message                         character*20
c            cinfo  : path of error                         character*80
c            info   : error number                               integer
c
c     CALL   :
c        GESTERR : gestion of errors
c
c-----------------------------------------------------------------------
c
      subroutine GESTOP ( iosort, finmes, cinfo, info )
c
      implicit none
c
      integer iosort, info
      character*20 finmes
      character*80 cinfo
c
c-----------------------------------------------------------------------
c
      if (info.ne.0) then
         call GESTERR ( cinfo, finmes )
         write(iosort,*)'**************error********************'
         write(iosort,*)'             error ',info,' in : ',cinfo
         write(iosort,*)'**************error********************'
         write(iosort,*)
         stop
      end if
c
      return
      end
c
c=======================================================================
c
c     subroutine VWIN
c
c     verification workspace - initialisation
c
c=======================================================================
c
c     INPUT :
c            ioecr  : canal de sortie                            integer
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c            siwmax : size max of integer workspace              integer
c            sdwmax : size max of double precision workspace     integer
c
c     OUTPUT :
c            iwork  : integer workspace vector(siwmax)           integer
c            dwork  : double precision workspace vector(sdwmax)   double
c
c-----------------------------------------------------------------------
c
      subroutine VWIN ( ioecr, siwork, sdwork, siwmax, sdwmax,
     &                  iwork, dwork )
c
      implicit none
c
      integer ioecr, siwork, sdwork, siwmax, sdwmax, iwveri
      double precision dwveri
c
      integer iwork(*)
      double precision dwork(*)
c
      parameter ( iwveri = 12345, dwveri = 123.45678 )
c
c-----------------------------------------------------------------------
c
      if ( (siwork+1) .gt. siwmax ) then
         write(ioecr,*)
     &             'ERROR: memoire de travail insuffisante (entiers) ',
     &             'iwork demande =',(siwork+1),'>',siwmax
         stop
      end if
      iwork(siwork+1) = iwveri
c
      if ( (sdwork+1) .gt. sdwmax ) then
         write(ioecr,*)
     &             'ERROR: memoire de travail insuffisante (doubles) ',
     &             'dwork demande =',(sdwork+1),'>',sdwmax
         stop
      end if
      dwork(sdwork+1) = dwveri
c
      return
      end
c
c=======================================================================
c
c     subroutine VWOUT
c
c     verification workspace - sortie
c
c=======================================================================
c
c     INPUT :
c            ioecr  : canal de sortie                            integer
c            siwork : size of integer workspace                  integer
c            sdwork : size of double precision workspace         integer
c            iwork  : integer workspace vector(siwmax)           integer
c            dwork  : double precision workspace vector(sdwmax)   double
c
c-----------------------------------------------------------------------
c
      subroutine VWOUT ( ioecr, siwork, sdwork, iwork, dwork )
c
      implicit none
c
      integer ioecr, siwork, sdwork, iwveri
      double precision dwveri, myzero
c
      integer iwork(*)
      double precision dwork(*)
c
      parameter ( iwveri = 12345, dwveri = 123.45678, myzero = 1.d-14 )
c
c-----------------------------------------------------------------------
c
      write(ioecr,*)
      write(ioecr,*)'taille iwork : ',siwork
      if ( iwork(siwork+1)-iwveri .eq. 0 ) then
         write(ioecr,*)'verification iwork OK'
      else
         write(ioecr,*)'!!! debordement iwork !!!'
      end if
c
      write(ioecr,*)'taille dwork : ',sdwork
      if ( abs(dwork(sdwork+1)-dwveri) .le. myzero ) then
         write(ioecr,*)'verification dwork OK'
      else
         write(ioecr,*)'!!! debordement dwork !!!'
      end if
c
      write(ioecr,*)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPETF
c
c     impression en-tete fonction
c
c=======================================================================
c
c     INPUT :
c            ioecr  : canal de sortie                            integer
c            nomfon : nom de la fonction                    character*10
c
c-----------------------------------------------------------------------
c
      subroutine IMPETF ( ioecr, nomfon )
c
      implicit none
c
      integer ioecr
      character*10 nomfon
c
c-----------------------------------------------------------------------
c
      if (ioecr.ne.6) then
         write(*,*)
         write(*,*)'_________________________________________________'
         write(*,*)
         write(*,*)'       ----- test de ',nomfon,' ----- '
         write(*,*)
      endif
c
      write(ioecr,*)
      write(ioecr,*)'_________________________________________________'
      write(ioecr,*)
      write(ioecr,*)'       ----- test de ',nomfon,' ----- '
      write(ioecr,*)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPETT
c
c     impression en-tete test
c
c=======================================================================
c
c     INPUT :
c            ioecr  : canal de sortie                            integer
c            nomfon : nom de la fonction de test            character*25
c
c-----------------------------------------------------------------------
c
      subroutine IMPETT ( ioecr, nomfon )
c
      implicit none
c
      integer ioecr
      character*25 nomfon
c
c-----------------------------------------------------------------------
c
      if (ioecr.ne.6) then
         write(*,*)
         write(*,*)'================================================='
         write(*,*)
         write(*,*)'=====   test de ',nomfon,'    ===== '
         write(*,*)
         write(*,*)'================================================='
         write(*,*)
      endif
c
      write(ioecr,*)
      write(ioecr,*)'================================================='
      write(ioecr,*)
      write(ioecr,*)'=====   test de ',nomfon,'   ===== '
      write(ioecr,*)
      write(ioecr,*)'================================================='
      write(ioecr,*)
c
      return
      end

