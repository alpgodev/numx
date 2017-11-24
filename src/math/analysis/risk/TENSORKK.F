c=======================================================================
c
c     subroutine TENSORKK.F                                   
c
c     Estimation of co-kurtosis[k] matrix, i.e. k-th tansor - matrix H[k]
c     
c                        n
c                    1  ---         _             _             _             _ 
c     H[k](i,j,l) = --- \   [X(t,i)-X(i)]*[X(t,j)-X(j)]*[X(t,l)-X(l)]*[X(t,k)-X(k)] 
c                    n  / 
c                       ---
c                       t=1
c     where
c     
c                  n
c       _      1  ---
c       X(j) = -  \   X(i,j) for j = 1,...,p
c              n  / 
c                 ---
c                 i=1
c
c-----------------------------------------------------------------------
      SUBROUTINE TENSORKK ( n, p, k, l, x, dwork, Hk, info)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n      : number of values                                integer
c       p      : number of asset(s)                              integer
c       k      : k-th H matix                                    integer
c       l      : l-th H matrix                                   integer
c       x      : values (n*p)                                     double
c
c     WORKSPACES
c       dwork : p                                                 double
c
c     OUTPUT 
c       Hk     : k-th tansor matrix (p*p)                         double
c       info   : diagnostic argument                             integer
c
c     CALL 
c       MCM
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, p, k, l, info
      DOUBLE PRECISION x(n, *), Hk(p,*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c
c     local variables
      INTEGER i, j, imat, pdrho
      DOUBLE PRECISION a, b, c, d, tmp, sum
c
c     intrinsic function
      INTRINSIC DFLOAT
c
c-----------------------------------------------------------------------
c
c     initializations 
      info = 0
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdrho = 1
c     pdrho : pointer for rho vector (p)              
c
c     Total size of dwork array = p
c
c-----------------------------------------------------------------------
c
c     test k value
      IF ((k .LT. 1) .OR. (k .GT. p)) THEN
        info = -1000
        RETURN
      ENDIF 
c     test l value
      IF ((l .LT. 1) .OR. (l .GT. p)) THEN
        info = -1000
        RETURN
      ENDIF 
c             
c     means return(s)
      CALL MCM ( n, p, x, dwork(pdrho) )
c      
      tmp = 1.0/DFLOAT(n)
      DO j = 1,p
         DO i = 1,j
            sum = 0.
            DO imat = 1,n
                a = x(imat,j) - dwork(pdrho + j - 1)
                b = x(imat,i) - dwork(pdrho + i - 1)
                c = x(imat,l) - dwork(pdrho + l - 1)
                d = x(imat,k) - dwork(pdrho + k - 1)
                sum = sum + a*b*c*d
             ENDDO
             Hk(i,j) = tmp * sum
             IF ( i .NE. j ) Hk(j,i) = Hk(i,j)
          ENDDO
      ENDDO            
c 
      RETURN
      END
