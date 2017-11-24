c=======================================================================
c
c     Fichier d'utilitaires                 
c
c     sous-programmes pour les operations sur les matrices d'entiers 
c     vectorisees, version turbo utilisant BLAS et LAPACK pour les 
c     routines precedees de # 
c
c=======================================================================
c
c     regles de nomination des sous-programmes :
c
c     1ere lettre : type d'operation
c                   A : arrangements
c                   C : conversion
c                   D : difference
c                   E : extremums
c                   I : initalisation
c                   K : kombien ?
c                   M : moyenne
c                   N : norme
c                   O : operations multiples
c                   P : produit
c                   R : transpose
c                   S : somme
c                   T : trace
c                   W : who ?
c                   X : produit scalaire
c                   Y : copie
c                   Z : symetrisation
c     lettres suivantes : operandes
c                   BOR: bornes
c                   CM : colonne d'une matrice
c                   D  : vecteur associe a une matrice diagonale
c                   I  : index
c                   EQ : equal
c                   EV : elements d'un vecteur
c                   GE : greater or equal
c                   GT : greater
c                   LE : lower or equal
c                   LM : ligne d'une matrice
c                   LT : lower
c                   M  : matrice (n*m)
c                   MAX: maximum
c                   MC : matrice carree (n*n)
c                   MCP: matrice carree partielle (n*n)
c                   MCT: matrice carree transposee (n*n)
c                   MD : matrice diagonale (n*n)
c                   MDS: matrice diagonale symetrique (n*(n+1)/2)
c                   MI : matrice d'entiers (n*m)
c                   MIN: minimum
c                   MP : matrice partielle (n*m)
c                   MS : matrice symetrique (n*(n+1)/2)
c                   MT : matrice transposee (m*n)
c                   OC : ordre croissant
c                   OD : ordre decroissant
c                   R  : fonction inverse (du meme nom sans"R")
c                   V  : vecteur (n)
c                   VI : vecteur d'entiers (n)
c                   VT : vecteur transpose (n)
c                   X  : scalaire
c
c        AVIOCI  : sort a vector of integers in increase order
c                  with index in old vector
c        AVIODI  : sort a vector of integers in decrease order
c                  with index in old vector
c        CMSIMCI : converting a vector of symmetric matrix of integers
c                  in vector of full square matrix of integers
c        CVIMDI  : converting a vector of integers
c                  in vectorized full diagonal matrix of integers
c        EMIBOR  : computing the minimum and maximum
c                  of a matrix of integers
c        EMIMAX  : computing the maximum of a matrix of integers
c        EMIMIN  : computing the minimum of a matrix of integers
c        EVIBOR  : computing the minimum and maximum
c                  of a vector of integers
c        EVIMAX  : computing the maximum of a vector of integers
c        EVIMIN  : computing the minimum of a vector of integers
c        IMI     : initialization at an integer of a matrix of integers
c        IMDI    : initialization at an integer of a diagonal matrix
c                  of integers
c        IVI     : initialization at an integer of a vector
c        KEVIEQ  : Computing how many elements of a vector
c                  of integer = value
c        SEVI    :  Sum of the element of vector.
c        WVIEQ   : who is equal to n in a vector ?
c                  gives a vector of indexes
c        WVIGE   : who is greater or equal than n in a vector ?
c                  gives a vector of indexes
c        WVIGT   : who is greater than n in a vector ?
c                  gives a vector of indexes
c        WVILE   : who is lower or equal than n in a vector ?
c                  gives a vector of indexes
c        WVILT   : who is lower than n in a vector ?
c                  gives a vector of indexes
c        YVI     : copy a vector of integers in a vector of integers
c        YVIP    : copy a part of a vector of integers
c                  in a vector of integers
c        YVXI  :   Initialization of a vector at an integer. 
c
c=======================================================================
c
c     subroutine AVIOCI
c
c     Sorts a vector of integers in increase order and returns indices
c
c-----------------------------------------------------------------------
      SUBROUTINE AVIOCI ( nvect, vecin, vecout, ind )
c-----------------------------------------------------------------------
c
c     INPUT
c       nvect  : vector size                                     integer
c       vecin  : input vector (nvect)                            integer
c
c     OUTPUT
c       vecout : sorted vector (nvect)                           integer
c       ind    : indices (nvect)                                 integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER nvect, ind(*), vecin(*), vecout(*)
      INTEGER i, j, indic, nout, a
c
c-----------------------------------------------------------------------
c
      ind(1) = 1
      vecout(1) = vecin(1)
      IF (nvect .GT. 1) THEN
         nout = 2
         indic = 1
         DO i=2,nvect
            a = vecin(i)
c
c           recherche du rang 
            indic = 1
            DO WHILE ( (indic.lt.nout) .and. (a.gt.vecout(indic)) )
               indic = indic + 1
            ENDDO
c
c           decalage
            DO j=nout,indic+1,-1
               vecout(j) = vecout(j-1)
               ind(j) = ind(j-1)
            ENDDO
            vecout(indic) = a
            ind(indic) = i
            nout = nout + 1
         ENDDO
      ENDIF
      RETURN
      END
c
c=======================================================================
c
c     subroutine AVIODI
c
c     Sorts a vector of integers (decrease order) and returns the 
c     indices in old vector
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vecin  : vector to sort (nvect)                     integer
c
c     OUTPUT 
c            vecout : sorted vector (nvect)                      integer
c            index  : indices of vector elements                 integer
c                     in the old vector (nvect)
c
c-----------------------------------------------------------------------
c
      subroutine AVIODI ( nvect, vecin, vecout, index )
c
      implicit none
c
      integer nvect
      integer index(*)
      integer vecin(*), vecout(*)
c
      integer i, j, indic, nout, a
c
c-----------------------------------------------------------------------
c
      index(1) = 1
      vecout(1) = vecin(1)
      if (nvect.gt.1) then
         nout = 2
         indic = 1
         do i=2,nvect
            a = vecin(i)
c
c           recherche du rang 
            indic = 1
            do while ( (indic.lt.nout) .and. (a.lt.vecout(indic)) )
               indic = indic + 1
            end do
c
c           decalage
            do j=nout,indic+1,-1
               vecout(j) = vecout(j-1)
               index(j) = index(j-1)
            end do
            vecout(indic) = a
            index(indic) = i
            nout = nout + 1
         end do
      endif
      return
      end
c
c=======================================================================
c
c     subroutine CMSIMCI
c
c     Converting a vector of symmetric matrix of integers
c     in vector of full square matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nmat   : size of square matrix mat                  integer
c            vects  : symmetric matrix as vect.(nmat*(nmat+1)/2) integer
c
c     OUTPUT 
c            vect   : matrix as vector(nmat*nmat)                integer
c
c-----------------------------------------------------------------------
c
      subroutine CMSIMCI ( nmat, vects, vect)
c
      implicit none
c
      integer nmat, vects(*), vect(nmat,*)
c
      integer i, j, k
c
c-----------------------------------------------------------------------
c
      k = 0
      DO j=1,nmat
         DO i=1,j
            k = k + 1
            vect(i,j) = vects(k)
            IF (i .NE. j) vect(j,i) = vects(k)
         ENDDO
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine CVIMDI
c
c     Converting a vector of integers
c     in vectorized full diagonal matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nmat   : size of matrix dmat and vector vect        integer
c            vect   : vector(nmat)                               integer
c
c     OUTPUT 
c            dmat   : diagonal matrix as vector(nmat*nmat)       integer
c
c-----------------------------------------------------------------------
c
      subroutine CVIMDI ( nmat, vect, dmat)
c
      implicit none
c
      integer nmat
      integer dmat(nmat,*), vect(*)
c
      integer i,j
c
c-----------------------------------------------------------------------
c
      do j=1,nmat
         do i=1,nmat
            dmat(i,j) = 0
         end do
         dmat(j,j) = vect(j)
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EMIBOR
c
c     Computing the minimum and maximum of a vect. matrix of integers
c
c=======================================================================
c
c     INPUT :
c            nrow   : number of rows of the matrix               integer
c            mcol   : number of columns of the matrix            integer
c            mat    : matrix(nrow*mcol)                          integer
c
c     OUTPUT :
c            mmin   : minimum of the matrix                      integer
c            mmax   : maximum of the matrix                      integer
c 
c-----------------------------------------------------------------------
c
      subroutine EMIBOR ( nrow, mcol, mat, mmin, mmax )
c
      implicit none
c
      integer nrow, mcol, mmin, mmax, mat(*)
c
      integer mcour, i
c
c-----------------------------------------------------------------------
c
      mmin = mat(1)
      mmax = mat(1)
      do i=2,nrow*mcol
         mcour = mat(i)
         if (mcour.lt.mmin) mmin = mcour
         if (mcour.gt.mmax) mmax = mcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EMIMAX
c
c     Computing the maximum of a vectorized matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nrow   : number of rows of the matrix               integer
c            mcol   : number of columns of the matrix            integer
c            mat    : matrix(nrow*mcol)                          integer
c
c     OUTPUT 
c            mmax   : maximum of the matrix                      integer
c 
c-----------------------------------------------------------------------
c
      subroutine EMIMAX ( nrow, mcol, mat, mmax )
c
      implicit none
c
      integer nrow, mcol, mmax
      integer mat(*)
c
      integer i, mcour
c
c-----------------------------------------------------------------------
c
      mmax = mat(1)
      do i=2,nrow*mcol
         mcour = mat(i)
         if (mcour.gt.mmax) mmax = mcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EMIMIN
c
c     Computing the minimum of a vectorized matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nrow   : number of rows of the matrix               integer
c            mcol   : number of columns of the matrix            integer
c            mat    : matrix(nrow*mcol)                          integer
c
c     OUTPUT 
c            mmin   : minimum of the matrix                      integer
c 
c-----------------------------------------------------------------------
c
      subroutine EMIMIN ( nrow, mcol, mat, mmin )
c
      implicit none
c
      integer nrow, mcol, mmin
      integer mat(*)
c
      integer i, mcour
c
c-----------------------------------------------------------------------
c
      mmin = mat(1)
      do i=2,nrow*mcol
         mcour = mat(i)
         if (mcour.lt.mmin) mmin = mcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EVIBOR
c
c     Computing the minimum and maximum of a vector of integers
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c
c     OUTPUT 
c            vmin   : minimum of the vector                      integer
c            vmax   : maximum of the vector                      integer
c
c-----------------------------------------------------------------------
c
      subroutine EVIBOR ( nvect, vect, vmin, vmax )
c
      implicit none
c
      integer nvect, vmin, vmax
      integer vect(*)
c
      integer i, vcour
c
c-----------------------------------------------------------------------
c
      vmin = vect(1)
      vmax = vect(1)
      do i=2,nvect
         vcour = vect(i)
         if (vcour.lt.vmin) vmin = vcour
         if (vcour.gt.vmax) vmax = vcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EVIMAX
c
c     Computing the maximum of a vector of integers
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c
c     OUTPUT 
c            vmax   : maximum of the vector                      integer
c
c-----------------------------------------------------------------------
c
      subroutine EVIMAX ( nvect, vect, vmax )
c
      implicit none
c
      integer nvect, vmax
      integer vect(*)
c
      integer i, vcour
c
c-----------------------------------------------------------------------
c
      vmax = vect(1)
      do i=2,nvect
         vcour = vect(i)
         if (vcour.gt.vmax) vmax = vcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine EVIMIN
c
c     Computing the minimum of a vector of integers
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c
c     OUTPUT 
c            vmin   : minimum of the vector                      integer
c
c-----------------------------------------------------------------------
c
      subroutine EVIMIN ( nvect, vect, vmin )
c
      implicit none
c
      integer nvect, vmin
      integer vect(*)
c
      integer i, vcour
c
c-----------------------------------------------------------------------
c
      vmin = vect(1)
      do i=2,nvect
         vcour = vect(i)
         if (vcour.lt.vmin) vmin = vcour
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMI
c
c     Initialization at an integer of a matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mat    : matrix to initialize vector(nlmat*ncmat)   integer
c            val    : value to initialize the matrix             integer
c
c     OUTPUT 
c            mat    : matrix initialized vector(nlmat*ncmat)     integer
c
c-----------------------------------------------------------------------
c
      subroutine IMI ( nlmat, ncmat, mat, val )
c
      implicit none
c
      integer nlmat, ncmat, val
      integer mat(nlmat,*)
c
      integer i,j
c
c-----------------------------------------------------------------------
c
      do j=1,ncmat
         do i=1,nlmat
             mat(i,j) = val
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMDI
c
c     Initialization at an integer of a diagonal matrix of integers
c
c=======================================================================
c
c     INPUT 
c            nmat   : size of diagonal matrix mat                integer
c            mat    : matrix to initialize vector(nmat*nmat)     integer
c            val    : value to initialize the matrix             integer
c
c     OUTPUT 
c            mat    : diagonal matrix    vector(nmat*nmat)       integer
c
c-----------------------------------------------------------------------
c
      subroutine IMDI ( nmat, mat, val )
c
      implicit none
c
      integer nmat, val
      integer mat(nmat,*)
c
      integer i,j
c
c-----------------------------------------------------------------------
c
      do j=1,nmat
         do i=1,nmat
             if (i.eq.j) then
                mat(i,i) = val
             else
                mat(i,j) = 0
             end if
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IVI
c
c     Initialization at an integer of a vector
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector to initialize (nvect)               integer
c            val    : value to initialize the vector             integer
c
c     OUTPUT 
c            vect   : vector initialized (nvect)                 integer
c
c-----------------------------------------------------------------------
c
      subroutine IVI ( nvect, vect, val )
c
      implicit none
c
      integer nvect, val, vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      do i=1,nvect
          vect(i) = val
      end do
      return
      end
c
c=======================================================================
c
c     subroutine KEVIEQ  : Computing how many elements of a vector
c                          of integer = value
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            sum    : sum of the elements of the vector          integer
c                     if vect(i) = value
c
c-----------------------------------------------------------------------
c
      subroutine KEVIEQ ( nvect, vect, val, sum )
c
      implicit none
c
      integer nvect, val, sum
      integer vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      sum = 0
      do i=1,nvect
         if ( vect(i).eq.val )  sum = sum + 1
      end do
      return
      end
c
c=======================================================================
c
c     subroutine SEVI: Sum of the elements of a vector.
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c
c     OUTPUT 
c            sum    : sum of the elements of the vector          integer
c                     
c-----------------------------------------------------------------------
c
      subroutine SEVI ( nvect, vect, sum )
c
      implicit none
c
      integer nvect, sum
      integer vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      sum = 0
      do i=1,nvect
         sum = sum + vect(i)
      end do
      return
      end
c
c=======================================================================
c
c     subroutine WVIEQ
c
c     who is equal at n in a vector ?
c     gives a vector of indexes
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements = value, vector(nind)    integer
c
c-----------------------------------------------------------------------
c
      subroutine WVIEQ ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind, val
      integer index(*), vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
         if ( vect(i).eq.val ) then
            nind = nind + 1
            index(nind) = i
         endif
      end do
      return
      end
c
c=======================================================================
c
c     subroutine WVIGE
c
c     who is greater or equal than n in a vector ?
c     gives a vector of indexes
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements >= value, vector(nind)   integer
c
c-----------------------------------------------------------------------
c
      subroutine WVIGE ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind, val
      integer index(*), vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
         if ( vect(i).ge.val ) then
            nind = nind + 1
            index(nind) = i
         endif
      end do
      return
      end
c
c=======================================================================
c
c     subroutine WVIGT
c
c     who is greater than n in a vector ?
c     gives a vector of indexes
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements > value, vector(nind)    integer
c
c-----------------------------------------------------------------------
c
      subroutine WVIGT ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind, val
      integer index(*), vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
         if ( vect(i).gt.val ) then
            nind = nind + 1
            index(nind) = i
         endif
      end do    
      return
      end
c
c=======================================================================
c
c     subroutine WVILE
c
c     who is lower or equal than n in a vector ?
c     gives a vector of indexes
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements =< value, vector(nind)   integer
c
c-----------------------------------------------------------------------
c
      subroutine WVILE ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind, val
      integer index(*), vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
         if ( vect(i).le.val ) then
            nind = nind + 1
            index(nind) = i
         endif
      end do   
      return
      end
c
c=======================================================================
c
c     subroutine WVILT
c
c     who is lower than n in a vector ?
c     gives a vector of indexes
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vect   : vector(nvect)                              integer
c            val    : value                                      integer
c
c     OUTPUT 
c            nind   : size of the indexes vector (max : nvect)   integer
c            index  : index of elements < value, vector(nind)    integer
c
c-----------------------------------------------------------------------
c
      subroutine WVILT ( nvect, vect, val, nind, index )
c
      implicit none
c
      integer nvect, nind, val
      integer index(*), vect(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      nind = 0
      do i=1,nvect
         if ( vect(i).lt.val ) then
            nind = nind + 1
            index(nind) = i
         endif
      end do    
      return
      end
c
c=======================================================================
c
c     subroutine YVI
c
c     Copy a vector of integers in a vector of integers
c
c=======================================================================
c
c     INPUT 
c            nvect  : size of the vector                         integer
c            vectin : input vector(nvect)                        integer
c
c     OUTPUT 
c            vecout : output vector(nvect)                       integer
c
c-----------------------------------------------------------------------
c
      subroutine YVI ( nvect, vectin, vecout )
c
      implicit none
c
      integer nvect
      integer vectin(*), vecout(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
      do i=1,nvect
          vecout(i) = vectin(i)
      end do
      return
      end
c
c=======================================================================
c
c     subroutine YVIP
c
c     Copy a part of a vector in a vector
c
c=======================================================================
c
c     INPUT 
c            nin    : size of the input vector vecin             integer
c            vecin  : input vector (nin)                         integer
c            rfirst : first row of vecin to copy                 integer
c            nout   : size of the output vector vecout           integer
c            vecout : initial output vector (nout)               integer
c
c     OUTPUT 
c            vecout : modified output vector (nout)              integer
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      subroutine YVIP ( nin, vecin, rfirst, nout, vecout, info )
c
      implicit none
c
      integer nin, nout, rfirst, info
      integer vecin(*), vecout(*)
c
      integer i
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      if ( rfirst+nout-1.le.nin ) then
         do i=1,nout
            vecout(i) = vecin(i+rfirst-1)
         end do
      else
         info = -2
      end if  
      return
      end
c      
c=======================================================================
c
c     subroutine YVXI
c
c     Initialization of a vector at an integer. 
c
c-----------------------------------------------------------------------
      SUBROUTINE YVXI ( n, vect, vali)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n    : size of the input vector vecin             integer
c            vect  : input vector (nin)                         integer
c            vali  :                   integer
c
c     OUTPUT 
c            vecout : modified output vector (nout)              integer
c            info   : = 0 successful exit                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER i, n, vect(*), vali
c
      DO i = 1,n
        vect(i) = vali
      ENDDO
      RETURN
      END
