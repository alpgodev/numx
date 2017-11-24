c=======================================================================
c
c     Fichier d'utilitaires                 
c
c     sous-programmes inpout/output pour les tests
c     matrice et vecteurs d'entiers
c
c        LECIV  : lecture d'un vecteur d'entiers
c        LECIM  : lecture d'une matrice vectorisee d'entiers
c        IMPIV  : imprime un vecteur d'entiers
c        IMPIM  : imprime une matrice vectorisee d'entiers
c        IMPIMC : imprime une matrice carree d'entiers (vectorisee)
c        IMPIMS : imprime une matrice vectorisee symetrique d'entiers
c        ARIMV  : converting a matrix of integers in vector of matrix
c
c=======================================================================
c
c     subroutine LECIV
c
c     Lecture d'un vecteur d'entiers
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de lecture                           integer
c            nvect  : taille du vecteur                          integer
c
c     OUTPUT :
c            vect   : vecteur (nvect)                            integer
c
c-----------------------------------------------------------------------
c
      subroutine LECIV ( iosort, nvect, vect )
c
      implicit none
c
      integer iosort, nvect
      integer vect(*)
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
c     subroutine LECIM
c
c     Lecture d'une matrice vectorisee d'entiers
c
c=======================================================================
c
c     INPUT :
c            iosort : canal de lecture                           integer
c            nlmat  : nombre de lignes de la matrice mat         integer
c            ncmat  : nombre de colonnes de la matrice mat       integer
c
c     OUTPUT :
c            mat    : matrice vectorisee (nlmat*ncmat)           integer
c
c-----------------------------------------------------------------------
c
      subroutine LECIM ( iosort, nlmat, ncmat, mat )
c
      implicit none
c
      integer iosort, nlmat, ncmat
      integer mat(nlmat,*)
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
c     subroutine IMPIV
c
c     Imprime un vecteur d'entiers
c
c=======================================================================
c
c     IMPIV ( iosort, texte, nvect, vect )
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nvect  : dimension du vecteur vect
c            vect   : vecteur de dimension nvect
c
c-----------------------------------------------------------------------
c
      subroutine IMPIV ( iosort, texte, nvect, vect )
c
      implicit none
c
      integer iosort, nvect
      integer vect(*)
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
    1 format(50i10)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPIM
c
c     Imprime une matrice d'entiers
c
c=======================================================================
c
c     INPUT  :
c            iosort : canal d'impression
c            texte  : commentaire
c            nlmat  : nombre de lignes de la matrice mat
c            ncmat  : nombre de colonnes de la matrice mat
c            mat    : matrice de dimension nlmat*ncmat
c
c-----------------------------------------------------------------------
c
      subroutine IMPIM ( iosort, texte, nlmat, ncmat, mat )
c
      implicit none
c
      integer iosort, nlmat, ncmat
      integer mat(nlmat,*)
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
    1 format(50i10)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPIMC
c
c     Imprime une matrice carree d'entiers (vectorisee)
c
c=======================================================================
c
c     IMPIM ( iosort, texte, nmat, mat )
c
c     INPUT  :
c            iosort : canal d'impression
c            texte  : commentaire
c            nmat   : taille de la matrice carree mat
c            mat    : matrice de dimension nlmat*ncmat
c
c-----------------------------------------------------------------------
c
      subroutine IMPIMC ( iosort, texte, nmat, mat )
c
      implicit none
c
      integer iosort, nmat
      integer mat(nmat,*)
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
    1 format(50i10)
c
      return
      end
c
c=======================================================================
c
c     subroutine IMPIMS
c
c     Imprime une matrice symetrique vectorisee d'entiers
c
c=======================================================================
c
c     INPUT :
c            iosort : canal d'impression
c            texte  : commentaire
c            nmat   : taille de la matrice mat
c            vmat   : matrice symetrique vectorisee (nmat*(nmat+1)/2)
c     CALL   :
c          CMSIMCI : converting a vector of symmetric matrix of integers
c                    in vector of full square matrix of integers
c
c-----------------------------------------------------------------------
c
      subroutine IMPIMS ( iosort, texte, nmat, vmat )
c
      implicit none
c
      integer iosort, nmat
      integer vmat(*)
      character*20 texte
c
      integer i, j
      integer mat(nmat,nmat)
c
c-----------------------------------------------------------------------
c
      call CMSIMCI ( nmat, vmat, mat)
c
      write(iosort,*)
      write(iosort,*) texte
      do i=1,nmat
         write(iosort,1)(mat(i,j),j=1,nmat)
      end do
      write(iosort,*)
c
    1 format(50i10)
c
      return
      end
c
c=======================================================================
c
c     subroutine ARIMV
c
c     Converting a matrix of integers in vector of matrix
c
c=======================================================================
c
c     INPUT :
c            nmax   : max of first dimension of the matrix mat   integer
c            nlmat  : number of raws in matrix mat               integer
c            ncmat  : number of columns in matrix mat            integer
c            mat    : matrix (nlmat*ncmat)                       integer
c
c     OUTPUT :
c            vect   : vector (nlmat*ncmat)                       integer
c
c-----------------------------------------------------------------------
c
      subroutine ARIMV ( nmax, nlmat, ncmat, mat, vect )
c
      implicit none
c
      integer nmax, nlmat, ncmat
      integer mat(nmax,*), vect(*)
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

