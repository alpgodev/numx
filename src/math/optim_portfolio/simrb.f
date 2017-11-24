c=======================================================================
c
c     SIMRB                              
c
c     This function implements the multi-volatility constraint simulator
c     cf. allocrb.f
c
c-----------------------------------------------------------------------
      SUBROUTINE simrb ( indic, simext, nvar, ydual, funct, grad,
     &                   iwork, dwork )
c-----------------------------------------------------------------------
c
c     INPUT 
c       indic  : = 4, to compute function and gradients          integer
c       simext : entry point of an external subroutine
c                provided by the user (not used here)
c       nvar   : number of variables (nb. class)                 integer
c       ydual  : dual solution lambda (nvar)                      double
c
c     OUTPUT 
c       funct  : dual function value                              double
c       grad   : gradient value                                   double
c
c     WORKSPACE 
c
c       iwork : 5*n + neq + 2*nin + 11                           integer
c        info  (diagnostic argument)     in iwork(1)
c        n     (number of assets)        in iwork(2) 
c        neq   (nb. equal. constraint)   in iwork(3)
c        nin   (nb. inequal. constraint) in iwork(4)
c        class (class definition)        in iwork(5,...,4+n) 
c
c       dwork : n*(7*n+20+2*neq+2*nin)+4*neq+6*nin+14+nvar        double
c        cov    (cov. matrix)
c        Q      (Q matrix)
c        kappa  (risk aversion coef.)
c        rho    (expected returns)
c        cinf   (lower bounds)
c        csup   (upper bounds)
c        sigma  (budget constraints)
c        C      (constraint matrix)
c        b      (constraint vector)    
c        ... + local dwork (QUAPRO + lagr)            
c     CALL  
c       YM, YMCPI, PMX, YMCPIR, qp, OVTMCV, XV, IMX
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
      EXTERNAL simext
c
c     i/o arguments
      INTEGER indic, nvar
      DOUBLE PRECISION funct, ydual(*), grad(*)
c
c     workspaces      
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables 
      INTEGER i, j, info, npk
      INTEGER n, neq, nin
      INTEGER pind, pic, piw, pdcov, pdcov1, pdcov2, pdcov3, 
     &        pdlin, pdrho, 
     &        pdcinf, pdcsup, pdwopt, pdlagr, pdsigma, pdc, pdb, pdw,
     &        pdQ, pdkappa
      DOUBLE PRECISION lambda, sigma, mean, var, vark, sum, COEF, ZERO
      PARAMETER (COEF = -0.5E0, ZERO = 0.E0)
c
c     external subroutines
      EXTERNAL YM, YMCPI, PMX, YMCPIR, qp, OVTMCV, XV, IMX, initfeas
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      n    = iwork(2)
      neq  = iwork(3)
      nin  = iwork(4)
c
c     pointers for integer work space  : iwork
c     ---------------------------------------- 
      pic  = 5
c     pic  : pointer for classes ( n )
      pind = pic + ( n )
c     pind : pointer for index vector (n)
      piw  = pind + ( n )
c     piw  : pointer for initfeas 
c               needs 3*n1 + 2*nin + neq + 1
c               so 3*n+2*nin+neq+6
c      
c     piw  : pointer for QP internal workspace, who needs
c            ( 3*n + 2*nin + neq + 1 )
c             with neqtot = neq
c                  nintot = nin + 1
c             
c     Total so size of iwork array = 5*n + neq + 2*nin + 11
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcov = 1
c     pdcov : pointer for cov. matrix (n*n)
c     so n*n
      pdQ   = pdcov + ( n*n)
c     pdQ   : pointer for Q matrix (n*n)
c     so 2*n*n
      pdkappa = pdQ + ( n*n )
c     pdkappa : pointer for risk aversion coef. (1)
c     so 2*n*n + 1 
      pdrho = pdkappa + ( 1 )
c     pdrho : pointer for rho vector (n)
c     so n*(2*n+1)+1
      pdcinf = pdrho + ( n )
c     pdcinf : pointer for lower bounds (n)
c     so n*(2*n+2)+1
      pdcsup = pdcinf + ( n )    
c     pdcsup : pointer for lower bounds (n)
c     so n*(2*n+3)+1
      pdsigma= pdcsup + ( n )
c     pdsigma: pointer for volatilities budgets (nvar)
c     so n*(2*n+2)+1+nvar
      pdc    = pdsigma + ( nvar )
c     pdc    : pointer for constraint matrix (n*(neq+nin))
c     so n*(2*n+2+neq+nin)+1+nvar
      pdb    = pdc + ( n*( neq + nin ))
c     pdb    : pointer for constraint vector (neq+nin)
c     so n*(2*n+2+neq+nin)+1+nvar+neq+nin
      pdcov1 = pdb + ( neq + nin )
c     pdcov1 : pointer for cov. matrix block (n*n)
c     so n*(3*n+2+neq+nin)+1+nvar+neq+nin
      pdcov2 = pdcov1 + ( n*n )
c     pdcov2 : pointer for cov. matrix block (n*n)
c     so n*(4*n+2+neq+nin)+1+nvar+neq+nin
      pdcov3 = pdcov2 + ( n*n )
c     pdcov2 : pointer for cov. matrix block (n*n)
c     so n*(5*n+2+neq+nin)+1+nvar+neq+nin
      pdlin  = pdcov3 + ( n*n )
c     pdlin  : pointer for linear part (n)
c     so n*(5*n+3+neq+nin)+1+nvar+neq+nin
      pdwopt = pdlin + ( n )
c     pdwopt : pointer for wopt (n)  
c     so n*(5*n+4+neq+nin)+1+nvar+neq+nin    
      pdlagr = pdwopt + ( n )
c     pdlagr : pointer for Lagrange mult. (n+neq+nin)
c     so n*(5*n+5+neq+nin)+1+nvar+2*neq+2*nin
      pdw    = pdlagr + ( n + neq + nin )
c     pdw    : pointer for initfeas 
c               needs (n+1)*(neq+nin+2*(n+1)+11)+neq+3*nin  
c     pdw    : pointer for QP internal workspace who needs
c              (n*n + 6*n + 2*nintot)
c               so pdw needs (n+1)*(neq+nin+2*(n+1)+11)+neq+3*nin  
c                       needs n*(neq+nin+2*n+15)+2*neq +4*nin+13
c
c
c     Total size of dwork array = n*n + n*n + 1 + n + n + n
c                               + nvar + n*( neq + nin ) + (neq + nin )
c                               + n*n + n*n + n*n + n + n 
c                               + (n + neq + nin) + n*(neq+nin+2*n+15)+2*neq +4*nin+13
c                               = n*(5*n+5+neq+nin)+1+nvar+2*neq+2*nin + n*(neq+nin+2*n+15)+2*neq +4*nin+13
c                               = n*(7*n+20+2*neq+2*nin)+4*neq+6*nin+14+nvar
c
c-----------------------------------------------------------------------
c
c     initialize cov1
      CALL IMX ( n, n, dwork(pdcov1), ZERO ) 
c
c     block risk budgets
      DO i = 1,nvar
         lambda = ydual(i)    
         CALL IVX ( n, iwork(pind), ZERO ) ! pind -> 0         
c
c        dimension of block i 
         npk = 0
         DO j = 1,n
            IF (iwork(pic + j - 1) .EQ. i) THEN         
               npk = npk + 1
               iwork(pind + npk - 1) = j
            ENDIF
         ENDDO
         IF (npk .GT. 1) THEN
c
            CALL IMX ( npk, npk, dwork(pdcov2), ZERO ) ! pdcov2 -> 0
c         
c           extract block i
            CALL YMCPI ( n, dwork(pdcov), npk, iwork(pind),
     &                   dwork(pdcov2), info )
c
c           product lambda(i)*cov(i)     
            CALL PMX (npk, npk, dwork(pdcov2), lambda, dwork(pdcov2))
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov2), n, iwork(pind),
     &                   dwork(pdcov1), info)
         ENDIF
      ENDDO
c            
c     product kappa*Q -> pdcov2
      CALL PMX (n, n, dwork(pdQ), dwork(pdkappa), dwork(pdcov2))
c
c     kappa*Q + Sum(Cov_i)      
      CALL SM ( n, n, dwork(pdcov1), dwork(pdcov2), dwork(pdcov3))
c
c     linear part: l(i) = -0.5*rho(i)
      CALL PVX ( n, dwork(pdrho), COEF, dwork(pdlin) )
c
c     Quadratic Solver
      CALL initfeas ( n, dwork(pdcinf), dwork(pdcsup), neq, nin,
     &                dwork(pdc), dwork(pdb), iwork(piw), dwork(pdw),
     &                dwork(pdwopt), info)
      IF (info .EQ. 1111) THEN
            funct = 1.e15
        RETURN
      ENDIF  
c     Initialization and feasibility check of the allocation problem
c     quadratic solver QP
      CALL qp ( n, dwork(pdcov3), dwork(pdlin), neq, nin,
     &          dwork(pdc), dwork(pdb), dwork(pdcinf), dwork(pdcsup),
     &          iwork(piw), dwork(pdw),
     &          dwork(pdlagr), dwork(pdwopt), info )
      
c     
c     objective function value
      CALL OVTMCV ( n, dwork(pdcov1), dwork(pdwopt), var)  ! w'*Cov*w 
      CALL OVTMCV ( n, dwork(pdcov2), dwork(pdwopt), vark) ! w'*Q*w
      CALL XV ( n, dwork(pdrho), dwork(pdwopt), mean)      ! rho'*w
c     (sigma**2)'*lambda
      sum = 0.0
      DO i = 1,nvar
        sigma = dwork(pdsigma + i - 1)
        sum = sum + sigma*sigma*ydual(i)
      ENDDO
      funct = mean + sum - var - vark
c      
c     gradient value (block risk budget): w'*Q(i)*w - sigma(i)
      DO i = 1,nvar
        sigma = dwork(pdsigma + i - 1)
c
c       dimension of block i
        npk = 0
        DO j = 1,n
            IF (iwork(pic + j - 1) .EQ. i) THEN
               npk = npk + 1
               iwork(pind + npk - 1) = j
            ENDIF
        ENDDO
        IF (npk .GT. 1) THEN
            CALL IMX ( npk, npk, dwork(pdcov2), ZERO ) ! cov2 <- 0
            CALL IMX ( n, n, dwork(pdcov1), ZERO )     ! cov1 <- 0
c           
c           extract block i
            CALL YMCPI ( n, dwork(pdcov), npk, iwork(pind),
     &                   dwork(pdcov2), info )
c
c           sub-block -> cov. matrix
            CALL YMCPIR (npk, dwork(pdcov2), n, iwork(pind),
     &                   dwork(pdcov1), info)
c
c           w'*Q(i)*w
            CALL OVTMCV ( n, dwork(pdcov1), dwork(pdwopt), var)
        ENDIF
        grad(i) = sigma*sigma - var
      ENDDO
      i=indic
      RETURN
      END
