c=======================================================================
c
c     subroutine CLUSTER                                     
c
c     Cluster Analysis 
c
c     Method: mixture approach based on multivariate Gaussian mixture 
c             model(s). The parameters are obtained from M.L. estimates
c             by the EM algorithm.
c
c-----------------------------------------------------------------------
      SUBROUTINE cluster ( n, d, k, x, iwork, dwork,
     &                     pk, nk, muk, covk, label, proba, f, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n        : number of points (n>1)                        integer
c       d        : dimension (d>0)                               integer
c       k        : number of class (k>0)                         integer
c       x        : data (n*d)                                     double
c
c     WORKSPACE
c       iwork    : d                                             integer
c       dwork    : d*(3*d + k + 9) + n                            double
c
c     OUTPUT 
c       pk       : proportion of each group (k)                   double
c       nk       : number of each group (k)                      integer
c       muk      : mean of each group (d*k)                       double
c       covk     : cov. matrix of each group (d)*(d*k)            double
c       label    : partition labels (n*k)                        integer
c       proba    : conditionale probability (n*k)                 double
c       f        : probability (n*k)                              double
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, d, k, nk(*), label(n,*), info
      DOUBLE PRECISION x(n,*), pk(*), proba(n,*), muk(d,*), covk(d,*), 
     &                 f(n,*)
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, l, itermax, piw, pdw, pdmuk, pdxi, ri
      DOUBLE PRECISION like, liketmp, variance, ZERO
      PARAMETER (ZERO = 0.0)
      REAL unif
c
c     external subroutines
      EXTERNAL em, YM
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      itermax = 10
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      piw = 1
c     piw : pointer for ESTEP, so ( d ) more
c
c     Total size of dwork array = ( d ) 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdw = 1
c     pdw : pointer for ESTEP, so ( d*(3*d + 9) ) more
      pdmuk = pdw + d*(3*d+9)
c     pdmuk : pointer for muk, so ( d*k ) more
      pdxi  = pdmuk + ( d*k )
c     pdxi  : pointer for x(i), so ( n ) more     
c
c     Total size of dwork array = d*(3*d + 9) 
c                               + d*k
c                               + n 
c                               = d*(3*d + k + 9) + n
c
c------------------------------------------------------------------------
c      open(unit=1,file='CLUSTER.txt',status='unknown')
c      write(1,*) "Begin CLUSTER"
c
c     parameter initialization 
c     ------------------------ 
c     mean initialization to x(1,k) and proportion to 1/K 
c
c     matrix identity by block (d)*(d*k)   
c
c            | 1 0 ... 0   ...   1 0 ... 0 |
c            | 0 1 ... 0   ...   0 1 ... 0 |
c     covk = | . .     .   ...   . .     . |
c            | . .     .   ...   . .     . |
c            | 0 0 ... 1   ...   0 0 ... 1 |
c
      DO l = 1,k
        pk(l) = 1./float(k)
        DO i = 1,d 
            CALL YCMV ( n, d, x, i, dwork(pdxi), info )
            CALL VARIAN ( n, dwork(pdxi), variance, info )
            covk(i, d*(l-1) + i) = variance
            muk(i,l) = x(l,i)
c            write(1,*) "mu initial(",i,",",l,")=", muk(i,l) 
        ENDDO    
      ENDDO    
c
c     initialisation by EM algorithm    
c      write(1,*) "Begin EM"
      
      CALL em ( itermax, n, d, k, x, iwork(piw), dwork(pdw),
     &          pk, nk, muk, covk, label, proba, f, like, info )   
      CALL YM ( d, k, muk, dwork(pdmuk))
c
c     initialization loop
      DO j = 1,100
c
c       initialisation: pk, covk, muk
        CALL IMX ( d, d*k, covk, ZERO )
        DO l = 1,k    
            pk(l) = 1./float(k)      
#ifdef INTELFOR 
        CALL RANDOM_NUMBER(HARVEST = unif)
#else
        unif = rand()
#endif            
            ri = int(unif*n)+1
            DO i = 1,d
                CALL YCMV ( n, d, x, i, dwork(pdxi), info )
                CALL VARIAN ( n, dwork(pdxi), variance, info )
                covk(i, d*(l-1) + i) = variance
                muk(i,l) = x(ri,i)
c                write(1,*) "mu Random=", muk(i,l)
            ENDDO
        ENDDO
c
c       EM algorithm
        CALL em ( itermax, n, d, k, x, iwork(piw), dwork(pdw),
     &            pk, nk, muk, covk, label, proba, f, liketmp, info )
c
        IF (liketmp .GT. like) THEN
            like = liketmp
            CALL YM ( d, k, muk, dwork(pdmuk))
c            write(1,*) "tests OK"
        ENDIF
     
      ENDDO
      
      CALL YM ( d, k, dwork(pdmuk), muk)
c
c     initialisation: pk (proportion), covk (cov. matrix)
      CALL IMX ( d, d*k, covk, zero )
      DO l = 1,k
        DO i = 1,d 
            covk(i, d*(l-1) + i) = 1.0
        ENDDO    
        pk(l) = 1./float(k)
      ENDDO
      itermax = 50
      CALL em ( itermax, n, d, k, x, iwork(piw), dwork(pdw),
     &          pk, nk, muk, covk, label, proba, f, like, info )
c
      RETURN
      END
     
