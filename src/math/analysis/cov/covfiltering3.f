c=======================================================================
c
c     subroutine COVFILRERING3
c
c     This function implements a cov. matrix correction for a filtering
c     level between 0% and 100%.
c
c-----------------------------------------------------------------------
      SUBROUTINE covfiltering3 ( p, cov, filter, iwork, dwork,
     &                           covopt, epsilon, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       p       : matrix size (p>0)                              integer
c       cov     : covariance matrix (p*p)                         double
c       filter  : filtering level (0<=filter<=1)                  double
c
c     WORKSPACE 
c       iwork   : 12*p + 3                                       integer 
c       dwork   : (p*(15*p + 61))/2 + 1                           double
c
c     OUTPUT 
c       covopt  : optimal cov. matrix (p*p)                       double
c       epsilon : calibration parameter                           double
c       info    : diagnostic argument                            integer
c
c     CALL   
c       SDLS    : Semidefinite Least Square optimization
c     
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER p, info
      DOUBLE PRECISION filter, epsilon, cov(*), covopt(*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER alpha, N, nbCst, pisdls, pdsdls, pdb, pdC, pdS
      DOUBLE PRECISION Lmin, Lmax, a, b, x
      DOUBLE PRECISION traceCov, cstPrecision, EPS, minEpsilon
      PARAMETER ( EPS = 1.E-15, minEpsilon = 1.E-8 ) 
      PARAMETER ( cstPrecision  = 1.E-15 )
c
c     external subroutines
      EXTERNAL TM, sdls, EVMIN, EVMAX, eigenvd
c
c     intrinsic functions
      INTRINSIC MAX
c
c-----------------------------------------------------------------------
c
c     initializations 
      info  = 0     ! diagnostic argument
      N     = 100   ! nb. filter increment
      nbCst = 0     ! number of SDLS constraints
      epsilon = 0.0 ! correction parameter
      CALL YM ( p, p, cov, covopt ) ! copy cov -> covopt    
c
c     pointers for integer work space : iworka
c     ---------------------------------------
      pisdls = 1
c     pisdls : pointer for SDLS (12*p + 3)
c
c     Total size of iwork array = 12*p + 3
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------   
      pdS    = 1
c     pdS    : pointer for eigenvalues (p)
      pdC    = pdS + ( p )
c     pdC    : pointer for cst. matrix SDLS (p*p)     
      pdb    = pdC + ( p*p )
c     pdb    : pointer for cst. vector SDLS (1)   
      pdsdls = pdb + ( 1 ) 
c     pdsdls : pointer for SDLS           
c                    5*p*(p+1)/2 + p*(4*p + 27)
c
c     Total size of dwork array = p + p*p + 1 + 5*p*(p+1)/2  + p*(4*p + 27)
c                               = (p*(15*p + 61))/2 + 1
c
c-----------------------------------------------------------------------
c
c     if filter=0 no correction
      IF (filter .LT. EPS) RETURN
c
c     cov. matrix trace
      CALL TM ( p, cov, traceCov )
      IF (traceCov .LT. EPS) RETURN
c
c     Singular value decomposition  
c     cov = U*Diag(S)*U', with sorted singular values S(i) >= S(i+1)
      CALL eigenvd (p, cov, dwork(pdsdls), dwork(pdS), 
     &              dwork(pdsdls + p*p + 10*p), info) 
c
c     min/max eigenvalue
      CALL EVMIN ( p, dwork(pdS), Lmin )
      CALL EVMAX ( p, dwork(pdS), Lmax )
c
c     SDP constraint (min eigenvalue)
      Lmin = MAX(minEpsilon, Lmin)
c
c     case Lmax = Lmin
      IF ((Lmax-Lmin) .LT. EPS) THEN
        info    = 1310 
        epsilon = Lmax
        RETURN
      ENDIF
c
c     exponent coefficient
      alpha = 6
c
c     a and b coefficients 
      a = (Lmax - Lmin)/(Lmax**alpha - Lmin**alpha)
      b = Lmin - a*(Lmin**alpha)
c
c     filter [0, 1] -> SDLS epsilon
      x = Lmin + 100*filter*(Lmax-Lmin)/N
      epsilon = a*(x**alpha) + b
c
c     optimization with SDLS
      CALL sdls ( p, cov, nbCst, dwork(pdC), dwork(pdb),
     &            cstPrecision, epsilon,
     &            iwork(pisdls), dwork(pdsdls), 
     &            covopt, info )
      IF (info .LT. 0) RETURN
c
      RETURN
      END
