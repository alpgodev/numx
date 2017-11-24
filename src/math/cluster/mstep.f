c=======================================================================
c
c     subroutine MSTEP                                     
c
c     M Step - Cluster Analysis 
c
c-----------------------------------------------------------------------
      SUBROUTINE mstep ( n, d, k, x, label, proba,
     &                   pk, nk, muk, covk, info )
c-----------------------------------------------------------------------
c
c     INPUT 
c       n        : number of points (n>1)                        integer
c       d        : dimension (d>0)                               integer
c       k        : number of class (k>1)                         integer
c       x        : data (n*d)                                     double
c       label    : partition label (n*k)                         integer
c       proba    : posteriori probability (n*k)                   double
c
c     OUTPUT 
c       pk       : proportion of each group (k)                   double
c       nk       : number of each group (k)                      integer
c       muk      : mean of each group (d*k)                       double
c       covk     : cov. matrix of each group (d)*(d*k)            double
c       info     : diagnostic argument                           integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, d, k, nk(*), label(n,*), info
      DOUBLE PRECISION x(n,*), pk(*), muk(d,*), proba(n,*), covk(d,*)
c
c     local variables
      INTEGER i, j, l, imat
      DOUBLE PRECISION EPS, sumd, invnk, a, b, p
      PARAMETER (EPS = 1.E-30)
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
c
c     M step (General Model), Sigma = lambda(k)*A'(k)D(k)A(k)
c     ---------------------
c
c     nk(k) : number of element(s) of the k-th group
      DO l = 1,k
        nk(l) = 0      
        DO i = 1,n
            nk(l) = nk(l) + label(i,l) 
        ENDDO
      ENDDO
c      
c     case nk(l)=0
      DO l = 1,k
        IF (nk(l) .EQ. 0) THEN 
             info = 6000
             RETURN
        ENDIF
      ENDDO  
c
c     pk(k) : proportion of the k-th group
      DO l = 1,k
        pk(l) = float(nk(l))/float(n)
      ENDDO
c
c     muk(d,k) : centre of the k-th group
      DO l = 1,k
        DO j = 1,d
            sumd = 0.0
            DO i = 1,n
                sumd = sumd + label(i,l)*x(i,j)
            ENDDO
            muk(j,l) = sumd/float(nk(l))
        ENDDO     
      ENDDO
c
c     covk(d,d*k) : cov matrix of the k-th group
      DO l = 1,k
c      
c       cov(d,d) : cov matrix of the k-th group     
        IF (nk(l) .EQ. 1) THEN
        DO i = 1,d
            DO j = 1,d
                covk(i, d*(l-1) + j) = 0.0
            ENDDO
            covk(i, d*(l-1) + i) = 1.E-15
        ENDDO    
        ELSE
        invnk = 1.0/float(nk(l))
c        
        DO j = 1,d
            DO i = 1,j
                sumd = 0.0
                DO imat = 1,n
                    p = proba(imat,l)
                    IF (p .GT. EPS) THEN
                        a = x(imat,i) - muk(i,l)
                        b = x(imat,j) - muk(j,l)
                        sumd = sumd + p*a*b
                    ENDIF    
                ENDDO
                covk(i, d*(l-1) + j) = invnk * sumd
                IF ( i .NE. j ) covk(j, d*(l-1)+i) = covk(i, d*(l-1)+j)
            ENDDO
        ENDDO
        ENDIF
      ENDDO
      RETURN
      END
