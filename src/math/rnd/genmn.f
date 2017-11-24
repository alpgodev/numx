c=======================================================================
c
c     GENMN                                                  
c 
c     Multivariate Normal random deviate
c
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------- 
      SUBROUTINE genmn(parm, x, work)
c-----------------------------------------------------------------------      
c
c     INPUT :
c      PARM --> Parameters needed to generate multivariate normal
c               deviates (MEANV and Cholesky decomposition of
c               COVM). Set by a previous call to SETGMN.
c               1 : 1                - size of deviate, P
c               2 : P + 1            - mean vector
c               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
c                                       decomposition of cov matrix
c                                             DOUBLE PRECISION PARM(*)
c
c      X  <-- Vector deviate generated (P)                       double
c      WORK <--> scratch array (P)                               double 
c
c     Method
c     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
c     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
c     3) trans(A)E + MEANV ~ N(MEANV,COVM)
c
c======================================================================
c
c     array arguments
      DOUBLE PRECISION parm(*), work(*), x(*)
c     
c     local scalars
      DOUBLE PRECISION ae
      INTEGER i, icount, j, p
c     
c     external functions
      DOUBLE PRECISION snorm
      EXTERNAL snorm
c     
c     intrinsic functions
      INTRINSIC int
c     
c     executable statements
      p = int(parm(1))
c
c     generate P independent normal deviates - WORK ~ N(0,1)
      DO 10,i = 1,p
          work(i) = snorm()
   10 CONTINUE
      DO 30,i = 1,p
c
c     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
c      decomposition of the desired covariance matrix.
c          trans(A)(1,1) = PARM(P+2)
c          trans(A)(2,1) = PARM(P+3)
c          trans(A)(2,2) = PARM(P+2+P)
c          trans(A)(3,1) = PARM(P+4)
c          trans(A)(3,2) = PARM(P+3+P)
c          trans(A)(3,3) = PARM(P+2-1+2P)  ...
c
c     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
          icount = 0
          ae = 0.0
          DO 20,j = 1,i
              icount = icount + j - 1
              ae = ae + parm(i+ (j-1)*p-icount+p+1)*work(j)
   20     CONTINUE
          x(i) = ae + parm(i+1)
   30 CONTINUE
      RETURN
      END
