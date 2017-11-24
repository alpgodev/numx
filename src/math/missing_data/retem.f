c=======================================================================
c
c     subroutine RETEM                                        
c
c     Expectation-Maximization (EM) algorithm 
c     -> recover missing observations in a time series 
c
c     Return data structure 
c
c           n : number of points
c           p : number of assets
c
c       R(n,p) = | x(1,1), ..., x(1,p) | 
c                | x(2,1), ..., -1000  |
c                | x(3,1), ..., x(3,p) |
c                |   .   , ...,   .    |
c                |   .   , ...,   .    |
c                | x(n,1), ..., x(n,p) |
c
c     where x(2,p) <= -1000 denote a missing return
c
c-----------------------------------------------------------------------
      SUBROUTINE retem ( n, p, X, miss, iwork, dwork, Y, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c        n      : max. number of values (n > 1)                  integer
c        p      : number of assets (p > 0)                       integer
c        X      : returns values (n*p)                            double
c        miss   : missing value (<= -1)                           double
c
c     WORKSPACE 
c        iwork  : 2*p                                            integer
c        dwork  : p*(13*p + 2*n + 10)                             double
c     
c     OUTPUT 
c        Y      : EM-recoved returns values (n*p)                 double
c        info   : diagnostic argument                            integer
c
c     CALL   
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, info
      DOUBLE PRECISION miss, X(n,*), Y(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER t, i, j, k, l, m, nobs, mt, nt, mtt, h
      INTEGER MaxLoop, CountLoop, Convergence
      INTEGER pindMiss, pindObs
      INTEGER pdX, pdy,
     &        pdrho, pdrhoN, pdrhoMiss, pdrhoObs,
     &        pdcov, pdcovN, pdcovm, 
     &        pdyMiss, pdyObs, pdM, pdyM, pdv,
     &        pdcovMO, pdcovMM, pdcovOO, pdinvCovOO, pdjms,
     &        pdCC, pdM1, pdM2, pdM3, pdMt, pdC
      DOUBLE PRECISION Tolerance, Distance, tmp, EPS, ZERO
      PARAMETER (EPS = 1.E-15, ZERO = 0.0)
c
c     external subroutines
      EXTERNAL IVX, YM, MCM, covm
c
c-----------------------------------------------------------------------
c
c     initialization
      info = 0
      CALL YM ( n, p, X, Y ) ! copy X -> Y
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pindMiss = 1
c     pindMiss : pointer missing data index (p)
      pindObs  = pindMiss + ( p )
c     pindObs  : pointer observed data index (p)           
c
c     Total size of iwork array  = 2*p      
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdX    = 1
c     pdX    : pointer for availaible NAV (n*p)
      pdrho  = pdX + ( n*p )
c     pdrho  : pointer for mean vector (p)
      pdrhoN = pdrho + ( p )
c     pdrhoN : pointer for new mean vector (p)
      pdcov  = pdrhoN + ( p )
c     pdcov  : pointer for covariance matrix (p*p)
      pdcovN = pdcov + ( p*p )
c     pdcovN : pointer for new covariance matrix (p*p)
      pdcovm = pdcovN + ( p*p ) 
c     pdcovm : pointer for COVM (p*n + p)  
      pdy    = pdcovm + ( n*p + p)
c     pdy    : pointer (p)
      pdyMiss = pdy + ( p )
c     pdyMiss : pointer missing data y (p) 
      pdyObs  = pdyMiss + ( p )
c     pdyObs  : pointer observed data y (p)
      pdM     = pdyObs + ( p )
c     
      pdyM    = pdM + ( p )
c
      pdv     = pdyM + ( p )
c
      pdrhoMiss  = pdv + ( p )
c     pdrhoMiss  : pointer missing mean (p)
      pdrhoObs   = pdrhoMiss + ( p )
c     pdrhoObs   : pointer observed mean (p)
      pdcovMO    = pdrhoObs + ( p )
c     pdcovMO    : pointer missing/observed cov (p*p)
      pdcovMM    = pdcovMO + ( p*p )
c     pdcovMM    : pointer missing/missing cov (p*p)
      pdcovOO    = pdcovMM + ( p*p )
c     pdcovOO    : pointer observed/observed cov (p*p)
      pdinvCovOO = pdcovOO + ( p*p )  
c     pdinvCovOO : pointer observed/observed inverse cov (p*p)
      pdCC       = pdinvCovOO + ( p*p )
c
      pdM1       = pdCC + ( p*p )
c
      pdM2       = pdM1 + ( p*p ) 
c
      pdM3       = pdM2 + ( p*p )
c
      pdMt       = pdM3 + ( p*p )
c
      pdC        = pdMt + ( p*p )
c
      pdjms      = pdC + ( p*p )
c     pdjms      : pointer for JMS (p*p)      
c
c     Total size of dwork array  = p*(13*p + 2*n + 10)
c             
c-----------------------------------------------------------------------
c     
c     E-M initialization
      h  = 1       ! horizon (log-returns
      m  = 0       ! number of missing data (one by row)
      DO t = 1,n
        DO i = 1,p
            IF (X(t,i) .LE. miss) THEN
                m = m + 1
                GOTO 100
            ENDIF         
        ENDDO
  100 CONTINUE        
      ENDDO
      IF (m .EQ. 0) RETURN ! if no missing data -> exit
      nobs = n - m         ! number of available dat
      IF (nobs .EQ. 0) THEN
        info = -9003 ! no enough available data for initialization 
        RETURN
      ENDIF
c     
      k = 0
      DO t = 1,n
        DO i = 1,p
            IF (X(t,i) .LE. miss) GOTO 101
        ENDDO
        k = k + 1
        DO i = 1,p
            dwork(pdX + nobs*(i-1) + k - 1) = X(t,i)
        ENDDO
  101 CONTINUE
      ENDDO
c
c     mean and cov. matrix
      CALL MCM (nobs,p,dwork(pdX),dwork(pdrho)) 
      CALL YV (p,dwork(pdrho),dwork(pdrhoN))                
      CALL covm (nobs,p,dwork(pdX),dwork(pdcovm),dwork(pdcov),info)
c     
      CALL YV (p*p,dwork(pdcov),dwork(pdcovN))
      IF (info .LT. 0) RETURN
c
c     EM loop
      tmp = dwork(pdrhoN) + SQRT(MAX(dwork(pdcovN),EPS))
      DO i = 2,p
        tmp = tmp + dwork(pdrhoN+i-1) + 
     &        SQRT(MAX(dwork(pdcovN+p*(i-1)+i-1),EPS))
      ENDDO
      Tolerance = 1.E-6*(tmp/(2.0*p)) ! algorithm tolerance
      MaxLoop     = 100  ! maximum EM loop
      CountLoop   = 0    ! loop counter
      Convergence = 0    ! convergence test (0 no, 1 yes)      
      DO WHILE ((Convergence .EQ. 0).AND.(CountLoop.LT.MaxLoop))
        CountLoop = CountLoop + 1  ! counter + 1
        CALL YV (p, dwork(pdrhoN), dwork(pdrho))
        CALL YV (p*p, dwork(pdcovN), dwork(pdcov))
c
c       Step 1: estimation
        CALL IVX (p*p, dwork(pdC), ZERO)
        CALL IVX (p*p, dwork(pdMt), ZERO) ! Mt = zeros(p,p)
        DO t = 1,n
             mt = 0 ! number of missing data at time t
             nt = 0 ! number of observed data at time t
             CALL IVX (p*p, dwork(pdCC), ZERO) ! CC = zeros(p,p)
             DO i = 1,p
                dwork(pdy + i - 1) = X(t,i)
                IF (X(t,i) .LE. miss) THEN
                    dwork(pdyMiss + mt)   = dwork(pdy+i-1)
                    iwork(pindMiss + mt)  = i ! index of missing data
                    dwork(pdrhoMiss + mt) = dwork(pdrho + i - 1)
                    mt = mt + 1 ! number of missing data at time t
                ELSE
                    dwork(pdyObs + nt)   = dwork(pdy+i-1) ! observation
                    iwork(pindObs + nt)  = i ! index of observed data
                    dwork(pdrhoObs + nt) = dwork(pdrho + i - 1)
                    nt = nt + 1 ! number of observed data at time t
                ENDIF
            ENDDO
c           
c           if all data are missing at date t -> exit             
            IF (mt .EQ. p) THEN
                info = -9001
                RETURN
            ENDIF     
            IF ((mt .GT. 0).AND. (nt .GT. 0)) THEN ! if missing data
c
c           construct cov. matrix for miss and observed data
c           cov[missing,obs] : vector of missing vs. observed covariances
c           cov[obs,obs]     : matrix of observed covariances
c           cov[missing,missing] : matrix of missing covariances
            DO i = 1,mt
                k = iwork(pindMiss + i - 1)
                DO j = 1,nt 
                    l = iwork(pindObs + j - 1) 
                    dwork(pdcovMO + mt*(j-1) + i-1) =
     &              dwork(pdcov + p*(k-1) + l-1)
                ENDDO
            ENDDO
            CALL YMCPI ( p, dwork(pdcov), mt, iwork(pindMiss),
     &                   dwork(pdcovMM), info)
            CALL YMCPI ( p, dwork(pdcov), nt, iwork(pindObs),
     &                   dwork(pdcovOO), info)
            CALL JMS (nt,dwork(pdcovOO),dwork(pdjms),dwork(pdinvCovOO), 
     &                info)
            IF (info .NE. 0) THEN
                info = -108 ! matrix is not sdp
                RETURN
            ENDIF
c
c           y(Miss) = M(Miss)+S(Miss,Obs)*inv(S(Obs,Obs))*(y(Obs)-M(Obs)
            CALL DV(nt,dwork(pdyObs),dwork(pdrhoObs),dwork(pdM))
            CALL PMV(nt,nt,dwork(pdinvCovOO),dwork(pdM),dwork(pdyM))
            CALL PMV(mt,nt,dwork(pdcovMO),dwork(pdyM),dwork(pdv))
            CALL SV(mt,dwork(pdrhoMiss),dwork(pdv),dwork(pdyMiss))
                        
c            
c           c(Miss,Miss) = S(Miss,Miss)-S(Miss,Obs)*inv(S(Obs,Obs))*S(Obs,Miss)
c           C(m,n) =  A(m,k) * B(k,n)
            CALL PM (mt,nt,nt,dwork(pdcovMO),dwork(pdinvCovOO),
     &                dwork(pdM1))
c           C(m,n) = A(m,k) * B'(n,k)     
            CALL PMMT (mt,nt,mt,dwork(pdM1),dwork(pdcovMO),dwork(pdM2))
c           C(m,n) = A(m,n) - B(m,n)   
            CALL DM (mt,mt,dwork(pdcovMM),dwork(pdM2),dwork(pdM3))
            DO i = 1,mt
                k = iwork(pindMiss + i - 1)
                DO j = 1,mt 
                    l = iwork(pindMiss + j - 1) 
                    dwork(pdCC + p*(k-1) + l - 1) =
     &              dwork(pdM3 + mt*(i-1) + j - 1) 
     
                ENDDO
            ENDDO
            CALL SV (p*p, dwork(pdCC), dwork(pdMt), dwork(pdC))
            CALL YV (p*p, dwork(pdC), dwork(pdMt))
c           
            mtt = 0
            DO i = 1,p
                IF (X(t,i) .LE. miss) THEN
                    mtt = mtt + 1
c                   replace missing data
                    dwork(pdy + mtt - 1) = 
     &              MAX(dwork(pdyMiss + mtt - 1),0.0)
	              Y(t,i) = dwork(pdy + mtt - 1)
	            ENDIF
      		ENDDO
         ENDIF
       ENDDO
c
c       Step 2: update
c       update mean rho_new = mean(Y)
        CALL MCM (n,p,Y,dwork(pdrhoN))
        tmp = 1.0/n 
c       update cov. matrix   
        CALL IVX (p*p, dwork(pdM1), ZERO) ! M1 = zeros(p,p)
        CALL IVX (p*p, dwork(pdM2), ZERO) ! M2 = zeros(p,p)
        CALL IVX (p*p, dwork(pdM3), ZERO) ! M3 = zeros(p,p)
        CALL PVX (p*p, dwork(pdC), tmp, dwork(pdM1))
        CALL COVM (n,p,Y,dwork(pdcovm),dwork(pdM2),info)
        tmp = (n-1.0)/n
        CALL PVX (p*p, dwork(pdM2), tmp, dwork(pdM3))
        CALL SV (p*p, dwork(pdM1), dwork(pdM3), dwork(pdcovN))
c    
c       test of algorithm convergence (mean and covariance matrix)
        CALL IVX (p, dwork(pdv), ZERO)    ! v = zeros(p)
        CALL IVX (p*p, dwork(pdM1), ZERO) ! M1 = zeros(p,p)
        CALL DV ( p, dwork(pdrho), dwork(pdrhoN), dwork(pdv) )
        CALL DV ( p*p, dwork(pdcov), dwork(pdcovN), dwork(pdM1) )
c       D4=[(M_new - M).^4; diag((S_new - S).^2) ]
        tmp = dwork(pdv) + dwork(pdM1)
        DO i = 2,p
            tmp = tmp + dwork(pdv+i-1)**4 + dwork(pdM1+p*(i-1)+i-1)**2
        ENDDO
c       Distance= mean(D4.^(1/4))
        Distance = SQRT(tmp/(2.0*p))
        Distance = SQRT(Distance)
        IF (Distance .LT. Tolerance) Convergence = 1
      ENDDO ! end while   
c      
      RETURN
      END
