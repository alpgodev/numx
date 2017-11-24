c=======================================================================
c
c     Random Generators (for Linux/Unix platform)            
c
c-----------------------------------------------------------------------
c
c        Fonctions List
c       
c        rn      : generate Normal(0,1)
c        rns     : generate Normal(mean,std)
c        rus     : generate Uniform(mean,half interval)
c
c     author: Yann Vernaz
c
c=======================================================================
c
c     function rn                                            
c
c     Standard Normal distribution generator N(0,1) 
c     Transformation method
c
c-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION rn()
c-----------------------------------------------------------------------
c     OUTPUT 
c       rn : generated normal deviate N(0,1)                      double
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     local variables
      DOUBLE PRECISION x1, x2, PI2, EPS
      PARAMETER ( PI2 = 6.283185307d0, EPS=1.d-15  )
      INTEGER init_seed_linux, s
c
c     external function      
      EXTERNAL init_seed_linux
c
c     intrinsic functions
      INTRINSIC SQRT, LOG, COS
c      
c     seeding function for UNIX FORTRAN random number generators
c     (cf. rne/random/init_seed_linux.f)
      s = init_seed_linux()           
c
      x1 = rand(s)
      DO WHILE (x1 .LT. EPS)
         x1 = rand()
      ENDDO
      x2 = rand()
      rn = SQRT(-2.0*LOG(x1))*COS(PI2*x2)
      RETURN
      END
c
c=======================================================================
c
c     function rns
c
c     Standard Normal distribution generator N(mean,std) 
c     with mean and standard deviation
c
c-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION rns ( mean, std )
c-----------------------------------------------------------------------      
c
c     INPUT 
c       mean : mean                                                double
c       std  : standard deviation                                  double
c
c     OUTPUT 
c       rns  : generated normal deviate N(mean,std)                double
c
c     CALL   
c       rn   : generate normal deviate N(0,1)     
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      DOUBLE PRECISION mean, std
c      
c     external functions
      DOUBLE PRECISION rn, snorm
      EXTERNAL rn, snorm
c
      rns = rn()*std + mean
c
      RETURN
      END
c
c=======================================================================
c
c     function rus                                           
c
c     Uniform distribution generator with mean and half interval
c
c-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION rus ( mean, hint )
c-----------------------------------------------------------------------      
c
c     INPUT 
c            mean   : mean of the normal distribution             double
c            hint   : half interval of the uniform distribution   double
c
c     OUTPUT 
c            rus    : scalar from an uniform distribution         double
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     i/o arguments
      DOUBLE PRECISION mean, hint
c
c     local variables      
      DOUBLE PRECISION x
      INTEGER init_seed_linux, s
c
c     external function      
      EXTERNAL init_seed_linux
c
c     seeding function for UNIX FORTRAN random number generators
c     (cf. rne/random/init_seed_linux.f)
      s = init_seed_linux()     
c
c     generate Uniform variate integer in the range [0,RAND_MAX] 
      x = rand(s)
c
c     Uniform variate with mean               
      rus = (x-0.5)*2.0*hint + mean
c
      RETURN
      END
