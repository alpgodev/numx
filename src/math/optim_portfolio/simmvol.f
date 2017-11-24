c=======================================================================
c
c     SIMMVOL                                 
c
c     This function implements the multi-volatility constraint simulator
c
c-----------------------------------------------------------------------
      SUBROUTINE simmvol ( indic, simext, nvar, ydual, funct, grad,
     &                     iwork, dwork )
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
c       iwork : 5*n + neq + 2*nin + 5                            integer
c        info  (diagnostic argument)     in iwork(1)
c        n     (number of assets)        in iwork(2) 
c        neq   (nb. equal. constraint)   in iwork(3)
c        nin   (nb. inequal. constraint) in iwork(4)
c        class (class definition)        in iwork(5,...,4+n) 
c
c       dwork : n*(4*n+neq+nin+12)+nvar+2*neq+4*nin              double
c        cov    (cov. matrix)
c        rho    (expected returns)
c        cinf   (lower bounds)
c        csup   (upper bounds)
c        sigma  (budget constraints)
c        C      (constraint matrix)
c        b      (constraint vector)        
c                     
c     CALL  
c       YM, YMCPI, PMX, YMCPIR, QP, OVTMCV, XV, IMX
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
      INTEGER i, j, l, k, info, n, neq, nin, npk
      INTEGER pind, pic, piw, pdcov, pdcov1, pdcov2, pdlin, pdrho, 
     &        pdcinf, pdcsup, pdwopt, pdlagr, pdsigma, pdc, pdb, pdw 
      DOUBLE PRECISION lambda, sigma, mean, var, sum, COEF, ZERO
      PARAMETER (COEF = -0.5E0, ZERO = 0.E0)
c
c     external subroutines
      EXTERNAL YM, YMCPI, PMX, YMCPIR, qp, OVTMCV, XV, IMX
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
c     piw  : pointer for QP (3*n + neq + 2*nin + 1)
c             
c     Total so size of iwork array = 5*n + neq + 2*nin + 5
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcov = 1
c     pdcov : pointer for cov. matrix (n*n)
      pdrho = pdcov + ( n*n )
c     pdrho : pointer for rho vector (n)
      pdcinf = pdrho + ( n )
c     pdcinf : pointer for lower bounds (n)
      pdcsup = pdcinf + ( n )    
c     pdcsup : pointer for lower bounds (n)
      pdsigma= pdcsup + ( n )
c     pdsigma: pointer for volatilities budgets (nvar)
      pdc    = pdsigma + ( nvar )
c     pdc    : pointer for constraint matrix (n*(neq+nin))
      pdb    = pdc + ( n*( neq + nin ))
c     pdb    : pointer for constraint vector (neq+nin)
      pdcov1 = pdb + ( neq + nin )
c     pdcov1 : pointer for cov. matrix block (n*n)
      pdcov2 = pdcov1 + ( n*n )
c     pdcov2 : pointer for cov. matrix block (n*n)
      pdlin  = pdcov2 + ( n*n )
c     pdlin  : pointer for linear part (n)
      pdwopt = pdlin + ( n )
c     pdwopt : pointer for wopt (n)      
      pdlagr = pdwopt + ( n )
c     pdlagr : pointer for Lagrange mult. (n+neq+nin)
      pdw    = pdlagr + ( n + neq + nin )
c     pdw    : pointer for QP (n*n + 6*n + 2*nin)
c
c     Total size of dwork array = n*(4*n+neq+nin+12)+nvar+2*neq+4*nin
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
c     linear part: l(i) = -0.5*rho(i)
      CALL PVX ( n, dwork(pdrho), COEF, dwork(pdlin) )
c
c     Quadratic Solver
      CALL qp ( n, dwork(pdcov1), dwork(pdlin), dwork(pdc), dwork(pdb),
     &             dwork(pdcinf), dwork(pdcsup), neq, nin, 
     &             iwork(piw), dwork(pdw),
     &             dwork(pdlagr), dwork(pdwopt), info )
     
c     
c     objective function value
      CALL OVTMCV ( n, dwork(pdcov1), dwork(pdwopt), var) ! w'*Q*w 
      CALL XV ( n, dwork(pdrho), dwork(pdwopt), mean )    ! rho'*w
c     (sigma**2)'*lambda
      sum = 0.0
      DO i = 1,nvar
        sigma = dwork(pdsigma + i - 1)
        sum = sum + sigma*sigma*ydual(i)
      ENDDO
      funct = mean + sum - var
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
      RETURN
      END
