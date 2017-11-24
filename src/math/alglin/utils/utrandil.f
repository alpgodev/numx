c=======================================================================
c
c     Integer Random Utilities                               
c
c     sous-programmes de generateur aleatoires et generateurs de bruit
c     sur des matrices ou vecteurs d'entiers
c
c-----------------------------------------------------------------------
c
c        irns    : Generates an integer single random deviate from a
c                  normal distribution with mean and standard deviation
c
c        irus    : Generates an integer single random deviate from an
c                  uniform distribution with mean and half interval
c
c        idie    : Generates a cast of die (n faces)
c
c        IDIEV   : Generates a vector from a cast of die (n faces)
c
c        IVIRN   : Initializes a vector of integers at an integer
c                  random from a normal distribution
c
c        IVIRU   : Initializes a vector of integers at an integer
c                  random from an uniform distribution
c
c        IMIRN   : Initializes a matrix of integers at an integer
c                  random from an normal distribution
c
c        IMSIRN  : Initializes a symmetric matrix of integers at
c                  an integer random from an normal distribution
c
c        IMIRU   : Initializes a matrix of integers at an integer
c                  random from an uniform distribution
c
c        IMSIRU  : Initializes a symmetric matrix of integers at
c                  an integer random from an uniform distribution
c
c=======================================================================
c
c     function irns
c
c     Generates an integer single random deviate from a
c     normal distribution with mean and standard deviation
c
c=======================================================================
c
c     INPUT :
c            mean   : mean of the normal distribution             double
c            std    : st_deviation of the normal distribution     double
c
c     OUTPUT :
c            irns   : integer from a normal distribution         integer
c
c     CALL   :
c            irn    : Generates a single random deviate from a normal
c                     distribution with mean at zero
c                     and standard deviation at one              integer
c
c-----------------------------------------------------------------------
      function irns ( mean, std )
      implicit none
      integer irns
      double precision mean, std, rn
      irns = int ( rn()*std + mean )
      return
      end
c
c=======================================================================
c
c     function irus
c
c     Generates an integer single random deviate from a
c     normal distribution with mean and half interval
c
c=======================================================================
c
c     INPUT :
c            mean   : mean of the normal distribution             double
c            hint   : half interval of the uniform distribution   double
c
c     OUTPUT :
c            irus   : integer from an uniform distribution       integer
c
c-----------------------------------------------------------------------
c
      function irus ( mean, hint )
      implicit none
      integer irus
      double precision mean, hint
      irus = int ( (rand()-0.5)*2.*hint + mean )
      return
      end
c
c=======================================================================
c
c     function idie
c
c     Generates a cast of die (n faces)
c
c=======================================================================
c
c     INPUT :
c            n      : number of faces of the die                 integer
c
c     OUTPUT :
c            idie   : cast of die                                integer
c
c-----------------------------------------------------------------------
      function idie ( n )
      implicit none
      integer n, idie
      idie = 0
      do while(idie.eq.0)
         idie = int ( rand() * n ) + 1
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IDIEV
c
c     Generates a vector from a cast of die (n faces)
c
c=======================================================================
c
c     INPUT :
c            nvect  : size of the vector                         integer
c            n      : number of faces of the die                 integer
c
c     OUTPUT :
c            vect   : nvect cast of die, vector(nvect)           integer
c
c     CALL   :
c        idie    : Generates a cast of die (n faces)
c
c-----------------------------------------------------------------------
c
      SUBROUTINE IDIEV ( nvect, n, vect )
      implicit none
      integer nvect, n, idie, i
      integer vect(*)
      do i=1,nvect
         vect(i) = idie(n)
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IVIRN
c
c     Initializes a vector of integers at a scalar
c     random from a normal distribution
c
c=======================================================================
c
c     INPUT :
c            nvect  : size of the vector                         integer
c            mean   : mean of the normal distribution             double
c            std    : st_deviation of the normal distribution     double
c
c     OUTPUT :
c            vnorm  : normal distribution vector(nvect)          integer
c
c     CALL   :
c            irns   : Generates an integer single random deviate from a
c                     normal distribution with mean
c                     and standard deviation
c
c-----------------------------------------------------------------------
c
      SUBROUTINE IVIRN ( nvect, mean, std, vnorm )
c
      implicit none
      integer nvect
      double precision mean, std
      integer vnorm(*)
      integer irns, i
      do i=1,nvect
         vnorm(i) = irns ( mean, std )
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IVIRU
c
c     Initializes a vector of integers at a scalar
c     random from an uniform distribution
c
c=======================================================================
c
c     INPUT :
c            nvect  : size of the vector                         integer
c            mean   : mean of the uniform distribution            double
c            hint   : half interval of the uniform distribution   double
c
c     OUTPUT :
c            vunif  : random vector(nvect)                       integer
c
c     CALL   :
c            irus   : Generates an integer single random deviate from an
c                     uniform distribution with mean and half interval
c
c-----------------------------------------------------------------------
c
      SUBROUTINE IVIRU ( nvect, mean, hint, vunif )
      implicit none
      integer nvect
      double precision mean, hint
      integer vunif(*)
      integer i, irus
      do i=1,nvect
         vunif(i) = irus ( mean, hint )
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMIRN
c
c     Initializes a matrix of integers at an integer
c     random from an normal distribution
c
c=======================================================================
c
c     INPUT :
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mean   : mean of the uniform distribution            double
c            std    : st_deviation of the normal distribution     double
c
c     OUTPUT :
c            mat    : random matrix as vector (nlmat*ncmat)      integer
c
c     CALL   :
c            irns   : Generates an integer single random deviate from a
c                     normal distribution with mean
c                     and standard deviation
c
c-----------------------------------------------------------------------
      SUBROUTINE IMIRN ( nlmat, ncmat, mean, std, mat )
      implicit none
      integer nlmat, ncmat
      double precision mean, std
      integer mat(nlmat,*)
      integer irns, i,j
      do i=1,ncmat
         do j=1,nlmat
            mat(j,i) = irns ( mean, std )
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMSIRN
c
c     Initializes a symmetric matrix of integers at an integer
c     random from an normal distribution
c
c=======================================================================
c
c     INPUT :
c            nmat   : size of the symmetric matrix mat           integer
c            mean   : mean of the uniform distribution            double
c            std    : st_deviation of the normal distribution     double
c
c     OUTPUT :
c            mat    : random symmetric matrix
c                     vector(nmat*(nmat+1)/2)                    integer
c
c     CALL   :
c            irns   : Generates an integer single random deviate from a
c                     normal distribution with mean
c                     and standard deviation
c
c-----------------------------------------------------------------------
      SUBROUTINE IMSIRN ( nmat, mean, std, mat )
      implicit none
      integer nmat
      double precision mean, std
      integer mat(*)
      integer irns, i,j,k
      k = 0
      do j=1,nmat
         do i=1,j
            k = k + 1
            mat(k) = irns ( mean, std )
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMIRU
c
c     Initializes a matrix of integers at an integer
c     random from an uniform distribution
c
c=======================================================================
c
c     INPUT :
c            nlmat  : number of rows of the matrix mat           integer
c            ncmat  : number of columns of the matrix mat        integer
c            mean   : mean of the uniform distribution            double
c            hint   : half interval of the uniform distribution   double
c
c     OUTPUT :
c            mat    : random matrix as vector (nlmat*ncmat)      integer
c
c     CALL   :
c            irus   : Generates an integer single random deviate from an
c                     uniform distribution with mean and half interval
c
c-----------------------------------------------------------------------
      SUBROUTINE IMIRU ( nlmat, ncmat, mean, hint, mat )
      implicit none
      integer nlmat, ncmat
      double precision mean, hint
      integer mat(nlmat,*)
      integer irus, i,j
      do i=1,ncmat
         do j=1,nlmat
            mat(j,i) = irus ( mean, hint )
         end do
      end do
      return
      end
c
c=======================================================================
c
c     subroutine IMSIRU
c
c     Initializes a symmetric matrix of integers at an integer
c     random from an uniform distribution
c
c=======================================================================
c
c     INPUT :
c            nmat   : size of the symmetric matrix mat           integer
c            mean   : mean of the uniform distribution            double
c            hint   : half interval of the uniform distribution   double
c
c     OUTPUT :
c            mat    : random symmetric matrix
c                     vector(nmat*(nmat+1)/2)                    integer
c
c     CALL   :
c            irus   : Generates an integer single random deviate from an
c                     uniform distribution with mean and half interval
c
c-----------------------------------------------------------------------
      SUBROUTINE IMSIRU ( nmat, mean, hint, mat )
      implicit none
      integer nmat
      double precision mean, hint
      integer mat(*)
      integer irus, i,j,k
      k = 0
      do j=1,nmat
         do i=1,j
            k = k + 1
            mat(k) = irus ( mean, hint )
         end do
      end do
      return
      end

