c=======================================================================
c
c     Random Utilities Functions              
c
c-----------------------------------------------------------------------
c
c        Fonctions List
c       
c        RNV     : Generates an history with normal random
c        RUV     : Generates an history with uniform random
c        BVN     : Generates V(i) = mean(i) + std*N(0,1)
c        BMN     : Generates M(i,j) = mean(i,j) + std*N(0,1)
c        BCMN    : Generates M(i,j) = mean(i,j) + std(j)*N(0,1)
c        BMSN    : M(i,j) = mean(i,j) + std*N(0,1) with M'=M (size n*(n+1)/2)
c        IVXRN   : Initializes a vector V(i) = N(mean,std), i=1,...,n
c        IVXRNS  : Initializes a vector V(i) = a*e(i)/sum(e(i)) 
c                  with e(i)~N(mean,std)
c        IVXRU   : Initializes a vector V(i) = Uniform(mean,hint), i=1,...,n
c        IVXRUS  : Initializes a vector V(i) = a*e(i)/sum(e(i)) 
c                  with e(i)~U(mean,hint)
c        IMXRN   : Initializes a matrix X(i,j) = N(mean,std)
c        IMSXRN  : Initializes a matrix X(i,j) = N(mean,std) and X'= X (symmetric)
c        IMXRU   : Initializes a matrix X(i,j) = Uniform(mean,hint)
c        IMSXRU  : Initializes a matrix X(i,j) = Uniform(mean,hint) 
c                  and X'= X (symmetric)
c
c=======================================================================
c
c     subroutine RNV                                         
c
c     Generates an history with a normal random on a vector
c
c-----------------------------------------------------------------------
      SUBROUTINE RNV ( n, vect, nmes, std, mnorm )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n      : vector size size of the vector                  integer
c       vect   : vector(n)                                        double
c       nmes   : number of measures with normal noise            integer
c       std    : standard deviation (n)                           double
c
c     OUTPUT 
c       mnorm  : matrix of initial vector with normal noise
c                     vectorized matrix(nmes*n)                   double
c
c     CALL   
c       IVXRN  : Initializes a vector at a scalar
c                     random from a normal distribution
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, nmes
      DOUBLE PRECISION vect(*), std(*), mnorm(nmes,*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
         CALL IVXRN ( nmes, vect(i), std(i), mnorm(1,i) )
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine RUV
c
c     Generates an history with an uniform random on a vector
c
c-----------------------------------------------------------------------
      SUBROUTINE RUV ( n, vect, nmes, hint, munif )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n  : vector size                         integer
c       vect   : vector(nvect)                               double
c       nmes   : number of measures with normal noise       integer
c       hint   : half interval of the uniform distribution   double
c                     vector(nvect)                               double
c
c     OUTPUT 
c       munif  : matrix of initial vector with uniform noise
c                     vectorized matrix(nmes*nvect)               double
c
c     CALL   
c       IVXRU  : Initializes a vector at a scalar
c                     random from an uniform distribution
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, nmes
      DOUBLE PRECISION vect(*), hint(*), munif(nmes,*)
c
c     local variables
      INTEGER i
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
         CALL IVXRU ( nmes, vect(i), hint(i), munif(1,i) )
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine BVN
c
c     Generates a normal random: V(i) = mean(i) + std*N(0,1)
c
c-----------------------------------------------------------------------
      SUBROUTINE BVN ( n, mean, std, V )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : vector size                                       integer
c       mean : mean vector (n)                                    double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       V  : random vector (n)                                    double
c
c     CALL   
c       rn : generated normal deviate N(0,1)  
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION std, mean(*), V(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rn
c
c     external function      
      EXTERNAL rn
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
         V(i) = mean(i) + std*rn() 
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine BMN
c
c     Generates: M(i,j) = mean(i,j) + std*N(0,1)
c
c-----------------------------------------------------------------------
      SUBROUTINE BMN ( n, p, mean, std, M )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : number of rows                                    integer
c       p    : number of columns                                 integer
c       mean : mean matrix (n*p)                                  double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       M    : random matrix (n*p)    double
c
c     CALL   
c       rn   : generated normal deviate N(0,1)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, p
      DOUBLE PRECISION std, mean(n,*), M(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION rn
c
c     external function      
      EXTERNAL rn
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
         DO j = 1,p
            M(i,j) = mean(i,j) + std*rn()
         ENDDO
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine BCMN                                       
c
c     Generates M(i,j) = mean(i,j) + std(j)*N(0,1)
c
c-----------------------------------------------------------------------
      SUBROUTINE BCMN ( n, p, mean, std, M )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       p    : number of columns                                 integer
c       mean : mean matrix (n*p)                                  double
c       std  : standard deviation vector (p)                      double
c
c     OUTPUT 
c       M    : random matrix (n*p)                                double
c
c     CALL   
c       rn   : generated normal deviate N(0,1)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, p
      DOUBLE PRECISION std(*), mean(n,*), M(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION rn, vol
c
c     external function      
      EXTERNAL rn
c
c-----------------------------------------------------------------------
c
      DO i = 1,p
         vol = std(i)
         DO j = 1,n
            M(j,i) = mean(j,i) + vol*rn()
         ENDDO
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine BMSN                                        
c
c     Generates M(i,j) = mean(i,j) + std*N(0,1) with M'=M
c
c-----------------------------------------------------------------------
      SUBROUTINE BMSN ( n, mean, std, M )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : matrix size                                       integer
c       mean : mean vestor (n*(n+1)/2)                            double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       M    : random symmetric matrix (n*(n+1)/2)                double
c
c     CALL   
c       rn   : generated normal deviate N(0,1)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION std, mean(*), M(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rn
c
c     external function
      EXTERNAL rn
c
c-----------------------------------------------------------------------
c
      DO i = 1,n*(n+1)/2
        M(i) = mean(i) + std*rn()
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine IVXRN                                      
c
c     Initializes a vector V(i) = N(mean,std), i=1,...,n
c
c-----------------------------------------------------------------------
      SUBROUTINE IVXRN ( n, mean, std, V )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : vector size                                       integer
c       mean : mean                                               double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       V    : random vector (n)                                  double
c
c     CALL   
c       rns  : generated normal deviate N(mean,std)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, std, V(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rns
c
c     external function
      EXTERNAL rns  
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
          V(i) = rns( mean, std )
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine IVXRNS                                    
c
c     Initializes a vector V(i) = a*e(i)/sum(e(i)) with e(i)~N(0,1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE IVXRNS ( n, mean, std, value, V)
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n     : vector size                                     integer
c       mean  : mean                                             double
c       std   : standard deviation                               double
c       value : value of the vector sum                          double
c
c     OUTPUT 
c        V    : normal vector(n)                                 double
c
c     CALL   
c       rns   : generated normal deviate N(mean,std)
c       PVX2  : scales a vector by a constant: a*V(n)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, std, value, V(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rns, sum
c
c     external functions
      EXTERNAL rns
c
c-----------------------------------------------------------------------
c
      sum = 0.
      DO i = 1,n
         V(i) = rns(mean,std)
         sum = sum + V(i)
      ENDDO
      sum = value / sum
      CALL PVX2 (n, V, sum)
      RETURN
      END
c
c=======================================================================
c
c     subroutine IVXRU                                      
c
c     Initializes a vector V(i) = Uniform(mean,hint), i=1,...,n
c
c-----------------------------------------------------------------------
      SUBROUTINE IVXRU ( n, mean, hint, V )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : vector size                                       integer
c       mean : mean of the uniform distribution                   double
c       hint : half interval of the uniform distribution          double
c
c     OUTPUT 
c       V    : uniform random vector (n)                          double
c
c     CALL   
c       rus  : generates a single random Uniform deviate
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, hint, V(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rus
c
c     external function 
      EXTERNAL rus
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
         V(i) = rus( mean, hint )
      ENDDO
      RETURN
      END
c
c=======================================================================
c
c     subroutine IVXRUS                                     
c
c     Initializes a vector V(i) = a*e(i)/sum(e(i)) with e(i)~U(mean,hint)
c
c-----------------------------------------------------------------------
      SUBROUTINE IVXRUS ( n, mean, hint, value, V )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n     : vector size                                      integer
c       mean  : mean                                              double
c       hint  : half interval of the uniform distribution         double
c       value : value of the vector sum                           double
c
c     OUTPUT 
c       V     : random vector (n)                                 double
c
c     CALL   
c       rus   : generates a single random Uniform deviate
c       PVX2  : scales a vector by a constant: a*V(n)
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, hint, value, V(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rus, sum
c
c     external function
      EXTERNAL rus      
c
c-----------------------------------------------------------------------
c
      sum = 0.
      DO i = 1,n
         V(i) = rus(mean,hint)
         sum = sum + V(i)
      ENDDO
      sum = value / sum
      CALL PVX2(n, V, sum)
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMXRN                                       
c
c     Returns a n-by-p matrix containing pseudo-random values
c     drawn from a Gaussian distribution
c     
c     X(i,j) = Normal(mean,std)
c
c-----------------------------------------------------------------------
      SUBROUTINE IMXRN ( n, p, mean, std, X )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       n    : number of rows                                    integer
c       p    : number of columns                                 integer
c       mean : mean                                               double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       X    : Gaussian random matrix (n*p)                       double
c
c     CALL   
c       rns  : generates a single Gaussian random deviate
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, p
      DOUBLE PRECISION mean, std, X(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION rns
c
c     external function
      EXTERNAL rns     
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
        DO j = 1,p
            X(i,j) = rns(mean,std)
        ENDDO
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMSXRN                                      
c
c     Returns a symmetric matrix containing pseudo-random values
c     drawn from a Gaussian distribution
c
c     X(i,j) = N(mean,std) and X'= X (symmetric)
c
c-----------------------------------------------------------------------
      SUBROUTINE IMSXRN ( n, mean, std, X )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       mean : mean of the uniform distribution                   double
c       std  : standard deviation                                 double
c
c     OUTPUT 
c       X    : Gaussian random symmetric matrix (n*(n+1)/2)       double
c
c     CALL   
c       rns  : generates a single Gaussian random deviate
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, std, X(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rns
c
c     external function
      EXTERNAL rns      
c
c-----------------------------------------------------------------------
c
      DO i = 1,n*(n+1)/2
        X(i) = rns(mean,std)
      ENDDO
c
      RETURN
      END

c
c=======================================================================
c
c     subroutine IMXRU                                      
c
c     Returns an n-by-p matrix containing pseudo-random values
c     drawn from a uniform distribution
c
c     X(i,j) = Uniform(mean,hint)
c
c-----------------------------------------------------------------------
      SUBROUTINE IMXRU ( n, p, mean, hint, X )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : number of rows                                    integer
c       p    : number of columns                                 integer
c       mean : mean of the uniform distribution                   double
c       hint : half interval of the uniform distribution          double
c
c     OUTPUT 
c       X    : Uniform random matrix (n*p)                        double
c
c     CALL   
c       rus  : generates a single random Uniform deviate
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n, p
      DOUBLE PRECISION mean, hint, X(n,*)
c
c     local variables
      INTEGER i,j
      DOUBLE PRECISION rus
c
c     external function
      EXTERNAL rus     
c
c-----------------------------------------------------------------------
c
      DO i = 1,n
        DO j = 1,p
            X(i,j) = rus(mean,hint)
        ENDDO    
      ENDDO
c
      RETURN
      END
c
c=======================================================================
c
c     subroutine IMSXRU                                      
c
c     Returns an symetric matrix containing pseudo-random values
c     drawn from a uniform distribution
c
c     X(i,j) = Uniform(mean,hint) and X'= X (symmetric)
c
c-----------------------------------------------------------------------
      SUBROUTINE IMSXRU ( n, mean, hint, X )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n    : matrix size                                       integer
c       mean : mean of the Uniform distribution                   double
c       hint : half interval of the Uniform distribution          double
c
c     OUTPUT 
c       X    : random symmetric matrix (n*(n+1)/2)                double
c
c     CALL   
c       rus  : generates a single random Uniform deviate
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      INTEGER n
      DOUBLE PRECISION mean, hint, X(*)
c
c     local variables
      INTEGER i
      DOUBLE PRECISION rus
c
c     external function
      EXTERNAL rus      
c
c-----------------------------------------------------------------------
c
      DO i = 1,n*(n+1)/2
        X(i) = rus(mean,hint)
      ENDDO
c
      RETURN
      END
