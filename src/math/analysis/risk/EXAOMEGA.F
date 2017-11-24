c=======================================================================
c
c     Function EXAOMEGA                                     
c
c     Ex-ante Omega function at threshold r     
c
c     Omega(r) = I(2)/I(1)
c
c     where
c
c                 / r                       /+oo
c                 |                         | 
c          I(1) = |   F(t)dt   and   I(2) = |  [1-F(t)]dt
c                 |                         | 
c                 /-oo                      /r
c  
c      F(.) is the cdf (cumulative distribution function) 
c
c      Omega ratio computed by Monte-Carlo Method
c
c-----------------------------------------------------------------------
      SUBROUTINE EXAOMEGA( n, w, cov, rho, r, 
     &                     iwork, dwork, omega, info)
c-----------------------------------------------------------------------
c
c     INPUT
c        n      : number of asset(s) (n>0)                        double
c        w      : n-array, portfolio weight(s)                    double 
c        cov    : n*n-array, covariance matrix (n by n)           double
c        rho    : n-array expected return(s)                      double
c        r      : threshold return                                double
c
c     WORKSPACE 
c        iwork  : vector ( 13*n )                                integer
c        dwork  : vector ( n*(5*n + 29) )                         double
c
c     OUTPUT
c        omega  : omega at threshold r                            double
c         info  : error code                                     integer
c
c     CALL   
c        EXARET : ex-ante return
c        EXARVO : ex-ante volatility
c
c----------------------------------------------------------------------
c     
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION r, omega, w(*), cov(*), rho(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c     
c     local variables
c      INTEGER seed, i 
      INTEGER m
      DOUBLE PRECISION num, den, EPS
      PARAMETER (EPS = 1.E-15)
c
c     external subroutines
      EXTERNAL EXARET, EXARVO    
c
c     external functions
      DOUBLE PRECISION rns
      EXTERNAL rns       
c     
c     intrinsic functions
      INTRINSIC max      
c
c-----------------------------------------------------------------------
c 
c     initialization
      omega = 0.0
      info  = 0  
      num   = 0.0
      den   = 0.0
c
c     number of Monte Carlo simulation  
      m     = 100
c
c     seed initialization for random
c#ifdef INTELFOR 
c      CALL TIME(seed)
c      CALL RANDOM_SEED(seed)
c#else
c      seed = time()
c      CALL srand(seed)
c#endif           
c
c     ex-ante mean 
c      CALL EXARET ( n, rho, w, meanx, info)
c      
c     ex-ante volatility
c      CALL EXARVO ( n, w, gamma, epskat, iwork, dwork,
c     &              volx, info)
c
c     Monte-Carlo m-loop
c      DO i = 1,m
c
c       generate random variable N(m,sigma)
c        rndx = rns ( meanx, volx ) 
c
c       numerator
c        num = num + max(rndx - r, 0.0)
c
c       denominator
c        den = den + max(r - rndx, 0.0)
c      ENDDO
c
c     omega ratio
c      IF (den .GT. eps) THEN
c        omega = num / den
c      ELSE
c        omega = 0.0
c        info = -1
c      ENDIF
c      
c     write omega
c      open(unit=1,file='omega.txt',status='unknown')
c	write(1,*) i
c	close(unit=1)
c
      RETURN
      END
