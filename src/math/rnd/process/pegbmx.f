c=======================================================================
c
c     subroutine PEGBMX                                      
c
c     Simulation of a Multidimentional Geometric Brownian Motion
c
c               dSt = St*[mu*dt + sigma*dWt]
c
c     Method: Euler discretization of SDE
c
c                S(1)    = S
c                S(1+t)  = S(t)*[1 + mu*dt +  sigma*sqrt(dt)*N(0,1)]
c
c                with dt = T/(N-1) 
c
c-----------------------------------------------------------------------
      SUBROUTINE pegbmx ( p, T, N, S, mu, cov, iwork, dwork, y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p          : number of element(s) (p > 1)                integer
c       T          : time (days, months, years)                  integer
c       N          : number of steps to maturity                 integer
c       S[p]       : initial values, (p)                          double
c       mu[p]      : process drift (p)                            double
c       cov[p,p]   : covariance matrix (p*p)                      double
c
c     WORKSPACE 
c       iwork      : 12*p                                        integer
c       dwork      : p(5*p + 31)                                  double
c
c     OUTPUT 
c       y[N+1,p]   : generated step, matrix(N*p)                  double
c       info       : diagnostic argument                         integer 
c
c     CALL   
c        YVLM      : copy a vector in a row of a vectorized matrix
c        GBMSTEP   : multidimensional GBM: dSt = mu*dt + chol(gamma)*dWt
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, T, N, info
      DOUBLE PRECISION S(*), mu(*), cov(*), y(N,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER pdstep, pdw, pdmu, pdcov, pichol, pdchol, 
     &        i, k, ks
      DOUBLE PRECISION dt, sqrdt
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     dt and sqrt(dt)
      dt    = (DFLOAT(T)/N)
      sqrdt = SQRT(dt)
c      
c     pointers for integer workspace : iwork
c     --------------------------------------
      pichol = 1
c     pichol : pointer for CHOL dwork (12*p)      
c
c     Total size of iwork array = 12*p
c
c     pointers for double precision work space : dwork
c     ------------------------------------------------
      pdstep = 1
c     pdstep : pointer for steps on each factor (p)
      pdmu   = pdstep + ( p )
c     pdmu   : pointer for mu*dt (p)
      pdcov  = pdmu + ( p )
c     pdcov  : pointer for Choleski matrix (p*p)
      pdchol = pdcov + ( p*p )
c     pdchol : pointer for CHOL dwork p(4*p + 27)      
      pdw    = pdchol + ( p*( 4*p + 27 ))
c     pdw   : pointer for GBMSTEP workspaces (2*p)
c
c     Total size of dwork array =  p + p + p*p
c                               +  p(4*p + 27)
c                               + 2*p
c                               = p(5*p + 31)
c
c-----------------------------------------------------------------------
c
c     initialization of trails at initial price
      CALL YVLM ( N, p, y, S, 1, info )
c
c     process drift: mu*dt
      CALL PVX ( p, mu, dt, dwork(pdmu) )
c
c     Choleski decomposition of covariance matrix
c      CALL CHOL ( p, cov, dwork(pdchol), dwork(pdcov), info )
c      IF (info .LT. 0) RETURN
c
c     robust choleski decomposition of covariance matrix
      CALL rcho ( p, cov, iwork(pichol),
     &            dwork(pdchol), dwork(pdcov), info )
      IF (info .LT. 0) RETURN
      
c
c     GBM trails: S(t+1) = S(t)*[1 + mu*dt + chol(gamma)*dWt]
      DO i = 1,N-1
c
c        simulation step 
         CALL GBMSTEP ( p, dwork(pdmu), dwork(pdcov), sqrdt, dwork(pdw),
     &                  dwork(pdstep))
c
c        price(n+1) = price(n)*(step+1)
         ks = pdstep
         DO k = 1,p
            y(i + 1,k) = y(i,k) * ( dwork(ks) + 1 )
            ks = ks + 1
         ENDDO
      ENDDO
c
      RETURN
      END
