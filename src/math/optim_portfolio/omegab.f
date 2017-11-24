c=======================================================================
c
c     subroutine OMEGAB                                     
c
c     Computes the variable omegaB for the RiskBudgetIT (cf. specification)
c     as well as the nbc values of (omegab^T Gamma_i omegab)
c
c-----------------------------------------------------------------------
      SUBROUTINE omegab ( n, cov, covb, nbc, class,
     &                    iwork, dwork,
     &                    omegabench, riskbench, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n          : portfolio size                              integer
c       cov        : covariance matrix (n*n)                      double
c       covb       : covariance assets-index (n)                  double
c       nbc        : nb. risk budgeting constraints (>=1)        integer 
c       class      : block definition (n)                        integer
c
c     WORKSPACE 
c       iwork      : 12*n+nbc                                    integer 
c       dwork      : n*(4*n+28)                                   double
c                    
c     OUTPUT 
c       omegabench : variable omega benchmark (n)                 double
c       riskbench  : nbc-vector of omegab^T Gamma_i omegab        double
c       info       : diagnostic argument                         integer
c
c     CALL   
c        IVI, IMX, SCHURS, OMCDMCT, OVTMCV, PMV
c
c-----------------------------------------------------------------------   
c      
      IMPLICIT NONE 
c
c     arguments i/o
      INTEGER n, info, nbc, class(*)
      DOUBLE PRECISION cov(*), covb(*), omegabench(*), riskbench(*)
c
c     worksoaces      
      DOUBLE PRECISION dwork(*) 
      INTEGER iwork(*)
      
c
c     local variables
      INTEGER izero, indk, indj, i, j, k, l, classSize, sumindex,
     &        pdmat, pdcovb, pdeigval, pdvpmat, pdschur, pdtemp,
     &        pdmatinv, piclas, pischur
      DOUBLE PRECISION epsil, temp, zero 
c
c     external function
      EXTERNAL IVI, IMX, SCHURS, OMCDMCT, OVTMCV, PMV      
c
c     initialization 
      epsil = 1.E-15
      izero = 0
      zero = 0.
      
      pdmat = 1
c needs n*n 
c so n*n
      pdcovb = pdmat +n*n
c needs n 
c so n*n + n = n*(n+1)
      
      pdeigval = pdcovb + n
c needs n 
c so n*(n+2)
      pdvpmat = pdeigval + n
c needs n*n 
c      so n*(2*n+2)
      pdmatinv = pdvpmat + n*n 
c needs n*n 
c      so n*(3*n+2)
      pdschur = pdmatinv + n*n 
c needs n*(n+26)
c so n*(4*n+28)
      pdtemp = pdmatinv + n*n
c needs n*n but is independant from schurs
c           and n*n < n*(n+26)
c       so   n*(4*n+28)
            
      piclas = 1
c needs nbc
c      
      pischur = piclas + nbc
c needs 12*n 
c      so 12*n+nbc

c     computation of the size of each class
      CALL ivi(nbc, iwork(piclas), izero)
      DO i = 1,n
        iwork(piclas+class(i)-1)=iwork(piclas+class(i)-1) + 1
      ENDDO 
      sumindex = 0
      DO i = 1,nbc
        indj = 0
        classSize = iwork(piclas+i-1)
        CALL IMX(classSize, classSize, dwork(pdmat), zero)
c       initialization of the matrix(classSize,classSize) 
c       corresponding to the i-th class
c       
        DO j = 1,n
            IF (class(j) .EQ. i) THEN
                indj = indj+1
                indk = 0
                DO k = 1,n
                    IF (class(k) .EQ. i) THEN
                        indk = indk+1
                        dwork(pdmat+indk-1+classSize*(indj-1))
     &                   = cov(k+(j-1)*n)
                    ENDIF
                ENDDO               
                dwork(pdcovb+indj-1) = covb(j)
            ENDIF
        ENDDO
c
c       Matrix Schur decomposition of a symmetric matrix
        CALL SCHURS(classSize, dwork(pdmat), iwork(pischur), 
     &              dwork(pdschur), dwork(pdeigval), dwork(pdvpmat), 
     &              info )
        IF (info .LT. 0) RETURN
c
c   
        DO l = 1,classSize
            temp = dwork(pdeigval+l-1)
            IF ((temp .GT. epsil) .OR. (temp .LT. -epsil)) THEN 
                dwork(pdeigval+l-1) = 1./temp
            ENDIF
        ENDDO
c
c     computation of P'*inv(Diag(sector(i))*P       
        CALL OMCDMCT(classSize, dwork(pdvpmat), dwork(pdeigval),
     &               dwork(pdtemp), dwork(pdmatinv))
        CALL PMV(classSize, classSize, dwork(pdmatinv), 
     &           dwork(pdcovb), omegabench(sumindex+1))
        CALL OVTMCV(classSize, dwork(pdmat), omegabench(sumindex+1), 
     &              riskbench(i))
        sumindex =  sumindex + classSize
      ENDDO
      RETURN 
      END
