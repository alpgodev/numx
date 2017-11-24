c=======================================================================
c
c     subroutine COVFILRERING
c
c     This function implements a cov. matrix correction for a filtering
c     level between 0% and 100%.
c
c-----------------------------------------------------------------------
      SUBROUTINE covfiltering ( p, cov, filter, iwork, dwork,
     &                          covopt, epsilon, info )
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
      INTEGER i, nbCst, pisdls, pdsdls, pdb, pdC, pdS
      DOUBLE PRECISION traceCov, epsMean, epsMedian, epsi, epsj, vpi, 
     &                 vpj, vEmp, vMean, vMedian, vmin, v, coef, a, b, 
     &                 cstPrecision, sum, EPS, 
     &                 minEPS
      PARAMETER ( EPS = 1.E-15, minEPS = 1.E-8, cstPrecision  = 1.E-15 )
c
c     external subroutines
      EXTERNAL YV, TM, sdls, eigenvd, calepsmedian
c
c     intrinsic functions
      INTRINSIC MAX, ABS
c
c-----------------------------------------------------------------------
c
c     initializations 
      info  = 0     ! diagnostic argument
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
      IF (filter .LE. EPS) RETURN
c
c     cov. matrix trace
      CALL TM ( p, cov, traceCov )
      IF (traceCov .LT. EPS) RETURN
c
c     Singular value decomposition   
c     cov = U*Diag(S)*U', with sorted singular values S(i) >= S(i+1)
      CALL eigenvd ( p, cov, dwork(pdsdls), dwork(pdS),
     &               dwork(pdsdls + p*p + 10*p), info )
c
c     &              dwork(pdsdls + p*p + 6*p), info)
c
c     eps mean and % exlpain variance
      epsMean = traceCov / p
c      CALL calepsmean(p, cov, epsMean, info)
      sum = 0.0
      DO i = 1,p
        epsi = dwork(pdS + i - 1)
        IF (epsi .GT. epsMean) THEN
            sum = sum + epsi 
        ENDIF
      ENDDO
      vMean = sum / traceCov
c
c     eps median and % explain variance
      CALL calepsmedian ( p, cov, dwork(pdsdls), epsMedian, info )
      sum = 0.0
      DO i = 1,p
        epsi = dwork(pdS + i - 1)
        IF (epsi .GT. epsMedian) THEN
            sum = sum + epsi 
        ENDIF
      ENDDO
      vMedian = sum / traceCov
c
c     empiric % explain variance between [2, 20] 
c      IF ( p .LE. 20) THEN
        a = -0.09999/18.0
        b = 0.9999 - 2.0*a
        vEmp = a*p + b  !
c      ELSE
c        vEmp = 0.9
c      ENDIF
c
c     % minimum of explain variance 
c     vmin = Max(vEmp, vMean, vMedian)
      vmin = MAX(vMean, vMedian)
      vmin = MAX(vEmp, vmin)
c
c     % explain variance: filter [0, 1] -> v [1 vmin]
      v = (vmin - 1.0)*filter + 1.0
c
c     minimum eigenvalue (linear interpolation)
      epsi = dwork(pdS)
      epsj = epsi
      vpi  = epsi/traceCov 
      vpj  = vpi 
      IF (p .GT. 1) THEN
        DO i = 2,p
            epsj = dwork(pdS + i - 1)
            vpj = vpj + epsj/traceCov
            IF (vpj .LT. v) THEN
                epsi = epsj
                vpi  = vpj
            ELSE
                GOTO 1
            ENDIF
        ENDDO
      ENDIF  
    1 CONTINUE
      IF (ABS(vpi - vpj) .LT. EPS) THEN
        coef = 0.0
      ELSE  
        coef = (v - vpj)/(vpi - vpj)
      ENDIF  
      epsilon = epsj + coef*(epsi - epsj)
      epsilon = MAX(epsilon, minEPS)
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
