c=======================================================================
c
c     subroutine MOMENT.F                                    
c
c     Estimation of k-th moment
c
c     m(k) = E[X^k] for k=1,2,3,4
c
c     ex. 3-th moment:
c               N   p
c              --- ----        
c     m(3) =   \   \    w(i)*w(j)*w(k)*X(t,i)*X(t,j)*X(t,k) 
c              /   /
c              --- ----
c              t=1 i,j,k
c
c-----------------------------------------------------------------------
      SUBROUTINE MOMENT ( n, p, k, w, x, kmoment, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c            n      : number of values                           integer
c            p      : number of asset(s)                         integer
c            k      : k-th moment (k=1,2,3,4)                     integer
c            w      : weights (p)                                 double
c            x      : values (n*p)                                double
c
c     OUTPUT 
c            kmoment: k-th moment                                 double
c            info   : diagnostic argument                        integer
c
c     CALL 
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, k, info
      DOUBLE PRECISION kmoment
      DOUBLE PRECISION w(*), x(n, *)
c
c     local variables
      INTEGER i, j, m, l, ll
      DOUBLE PRECISION tmp, sum
c
c-----------------------------------------------------------------------
c
c     initialization
      kmoment = 0.0
      sum = 0.0
c
c     test k value
      IF ((k .lt. 1).or.(k.gt.4)) THEN
        info = -1000
        RETURN
      ENDIF        
c
c     1-th moment   
      IF (k .eq. 1) THEN
      sum = 0.0
      DO i = 1,n
        tmp = 0.0
        DO j = 1,p
           tmp = tmp + w(j)*x(i,j)
        ENDDO    
         sum = sum + tmp
      ENDDO
      ENDIF
c
c     2-th moment   
      IF (k .eq. 2) THEN
      sum = 0.0
      DO i = 1,n
        tmp = 0.0
        DO j = 1,p
            DO l = 1,p
                tmp = tmp + w(j)*w(l)*x(i,j)*x(i,l)
            ENDDO  
        ENDDO    
        sum = sum + tmp
      ENDDO
      ENDIF      
c
c     3-th moment      
      IF (k .eq. 3) THEN
      sum = 0.0
      DO i = 1,n
        tmp = 0.0
        DO j = 1,p
            DO m = 1,p
                DO l = 1,p
                    tmp = tmp + w(j)*w(m)*w(l)*x(i,j)*x(i,m)*x(i,l)
                ENDDO    
            ENDDO
         ENDDO
         sum = sum + tmp
      ENDDO
      ENDIF
c
c     4-th moment      
      IF (k .eq. 4) THEN
      sum = 0.0
      DO i = 1,n
        tmp = 0.0
        DO j = 1,p
            DO m = 1,p
                DO l = 1,p
                    DO ll = 1,p
                          tmp = tmp + w(j)*w(m)*w(l)*w(ll)*x(i,j)*x(i,m)
     &                                *x(i,l)*x(i,ll)
                    ENDDO 
                ENDDO    
            ENDDO
         ENDDO
         sum = sum + tmp
      ENDDO
      ENDIF      
c
      kmoment = sum/n
c 
      RETURN
      END
