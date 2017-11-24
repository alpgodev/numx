c=======================================================================
c
c     Function EXPOMEGA                                      
c
c     Ex-post Omega function at threshold r     
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
c      Omega(r) is computed by Monte-Carlo method
c
c-----------------------------------------------------------------------
      SUBROUTINE EXPOMEGA( n, x, r, omega, info)
c-----------------------------------------------------------------------
c
c     INPUT
c        n      : number of return(s) (n > 1)                    integer
c        x      : return(s) (n)                                   double
c        r      : threshold return                                double
c
c     OUTPUT
c        omega  : omega at threshold r                            double
c         info  : error code                                     integer
c
c     CALL   
c       EXPRET  : ex-post mean return
c       VOLAT   : ex-post volatility
c
c----------------------------------------------------------------------
c     
c     arguments i/o
      INTEGER n, info
      DOUBLE PRECISION x(*), r, omega 
c     
c     local variables
      INTEGER seed, i, m
      DOUBLE PRECISION rndx, meanx, volx, num, den, EPS
      PARAMETER (EPS = 1.E-30)
c
c     external subroutines
      EXTERNAL EXPRET, VOLAT  
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
      m     = 1000 
c
c     seed initialization for random
#ifdef INTELFOR 
      CALL TIME(seed)
      CALL RANDOM_SEED(seed)
#else
      seed = time()
      CALL srand(seed)
#endif        
c
c     ex-post mean return
      CALL EXPRET ( n, x, meanx, info )
      IF (info .LT. 0) RETURN
c
c     ex-post volatility
      CALL VOLAT ( n, x, volx, info )
      IF (info .LT. 0) RETURN          
c
c     Monte-Carlo m-loop
      DO i = 1,m
c
c       generate random variable N(m,sigma)
        rndx = rns ( meanx, volx )
c
c       numerator (ratio) 
        num = num + max(rndx - r, 0.0) 
c
c       denominator (ratio)       
        den = den + max(r - rndx, 0.0)
      ENDDO
c
c     omega ratio 
      IF (den .GT. EPS) THEN
        omega = num / den
      ELSE
        omega = 0.0
c        info = -1
      ENDIF
c      
c     write omeaga
c     open(unit=1,file='omega.txt',status='unknown')
c	write(1,*) omega
c	close(unit=1)
c
      RETURN
      END
