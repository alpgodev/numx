c=======================================================================
c
c     subroutine ESTEP                                     
c
c     E Step - Cluster Analysis 
c
c-----------------------------------------------------------------------
      SUBROUTINE estep ( n, d, k, x, pk, muk, covk,
     &                   iwork, dwork, label, proba, f, like, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n        : number of points (n>1)                        integer
c       d        : dimension (d>0)                               integer
c       k        : number of class (k>0)                         integer
c       x        : data (n*d)                                     double
c       pk       : proportion of each group (k)                   double
c       muk      : mean of each group (d*k)                       double
c       covk     : cov. matrix of each group (d)*(d*k)            double
c
c     WORKSPACE
c       iwork    : d                                             integer
c       dwork    : d*(3*d + 9)                                    double
c
c     OUTPUT 
c       label    : partition label (n*k)                         integer
c       proba    : conditionale probability (n*k)                 double
c       f        : probability value (n*k)                        double
c       like     : likelihood value                               double       
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, d, k, label(n,*), info
      DOUBLE PRECISION x(n,*), pk(*), muk(d,*), proba(n,*), covk(d,*),
     &                 f(n,*), like 
c
c     workspaces
      INTEGER iwork(*)
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, l, row, col, pdcov, pdmu, pdxi, pipdf, pdpdf
      DOUBLE PRECISION normpdf, sum, EPS
      PARAMETER (EPS = 1.E-30)
c
c     external subroutines
      EXTERNAL YMP, YLMV, EVMAXI, getMNormalPdf      
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     pointers for integer work space  : iwork
c     ----------------------------------------
      pipdf = 1
c     pipdf : pointer for getMNormPdf, so ( d ) more
c
c     Total size of dwork array = ( d ) 
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdcov = 1
c     pdcov : pointer for cov. matrix, so ( d*d ) more
      pdmu  = pdcov + ( d*d )
c     pdmu  : pointer for mean vector, so ( d ) more
      pdxi  = pdmu + ( d )
c     pdxi  : pointer for x(i) point, so ( d ) more
      pdpdf = pdxi + ( d )
c     pdpdf : pointer for getMNormPdf, so ( d*(2*d + 7) ) more
c
c     Total size of dwork array = d*(2*d + 7) 
c
c------------------------------------------------------------------------
c      open(unit=1,file='E-STEP.txt',status='unknown')
c      write(1,*) "Begin E-Step"
c
c     E step (expectation)
c     --------------------
c
      DO l = 1,k
c
c       copy a part of covk (d)*(n*d) in cov (d*d)
        row = 1
        col = 1 + d*(l-1)
c
        CALL YMP ( d, d*k, covk, d, d, row, col, dwork(pdcov), info ) 
c
c       copy muk(:,l) in mu(d)
        CALL YCMV ( d, k, muk, l, dwork(pdmu), info )   
             
        DO i = 1,n
c
           CALL YLMV ( n, d, x, i, dwork(pdxi), info )
c
c          Normal probability density function (pdf)
c
c          returns the normal pdf with mean MU, and covariance matrix COV, 
c          at the values in x (size 1 by d). normpdf is scalar 
           CALL getMNormalPdf(d,dwork(pdxi),dwork(pdmu),dwork(pdcov),
     &                        iwork(pipdf),dwork(pdpdf),
     &                        normpdf, info)
           f(i,l) = normpdf
c
        ENDDO     
      ENDDO
c
c     conditional probabilities and likelihood
      like = 0.0
      DO i = 1,n
c
c       compute sum[p(k)*f(x(i),k))]
        sum = 0.0
        DO l = 1,k
            sum = sum + pk(l)*f(i,l)
        ENDDO
        DO l = 1,k         
            proba(i,l) = pk(l)*f(i,l)/sum   
        ENDDO
        like = like + LOG(sum)  
      ENDDO 
c
c     labels  
      DO i = 1,n         
c       copy the i-th row of c(n,k)
        CALL YLMV ( n, k, proba, i, dwork(pdmu), info )
c        
c        CALL YLMV ( n, k, f, i, dwork(pdmu), info )
c       max. of i-th row                  
        CALL EVMAXI ( k, dwork(pdmu), sum, col )
        DO l = 1,k
            label(i,l) = 0
        ENDDO
        label(i,col) = 1
      ENDDO  
c
      RETURN
      END
     
