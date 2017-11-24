c=======================================================================
c
c     subroutine BUILDTHGEN
c
c-----------------------------------------------------------------------
      SUBROUTINE buildthgen ( n, qlin, rho, w0, Q, ndual, ydual,
     &                        dwork, cov, lin)
c-----------------------------------------------------------------------
c
c     INPUT 
c       n       : number of values                               integer
c       qlin    : vector(n) linear part of the objective funtcion double
c       rho     : mean returns (n)                                double 
c       w0      : vector(n) second order portfolio                double
c       Q       : semi definite matrix (n*n)                      double
c       ndual   : number of variables                            integer
c       ydual   : dual solution lambda (ndual)                    double 
c
c     OUTPUT 
c       cov     : semi definite matrix (n*n)                      double
c       lin     : vector(n) linear part of the objective funtcion double
c
c     WORKSPACE 
c
c       dwork   : 3*n                                             double
c           
c     CALL  
c       YM, YV, PVX, SV
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o      
      INTEGER n, ndual
      DOUBLE PRECISION qlin(*), rho(*), w0(*), Q(*), ydual(*), cov(*),
     &                 lin(*)
c
c     workspaces
      DOUBLE PRECISION dwork(*)
c     
c     local variables      
      INTEGER i, pdx, pdy, pdz
      DOUBLE PRECISION a
c
c-----------------------------------------------------------------------     
c
c     pointers for double precision work space  : dwork
c     -------------------------------------------------
      pdx = 1
c     needs n
      pdy = pdx + n
c     needs n
      pdz = pdy + n
c     needs n
c     
c     Total size of dwork array = 3*n
c 
c-----------------------------------------------------------------------
c
c     initialize cov
      CALL YM(n, n, Q, cov)     ! cov. matrix
      DO i = 1,n
        cov((i-1)*n+i) = cov((i-1)*n+i) + 2*ydual(2)
      ENDDO
c
c     cov = Q + 2*ydual(1)*Id_n
      a = -ydual(1) 
      CALL PVX(n, rho, a, dwork(pdx))          ! pdx = a*rho
      a = -2*ydual(2)
      CALL PVX(n, w0, a, dwork(pdy))           ! pdy = a*w0
      CALL SV(n, qlin, dwork(pdx), dwork(pdz)) ! pdz = qlin + pdx
      CALL SV(n, dwork(pdz), dwork(pdy), lin)  ! lin = pdz + pdy
      RETURN
      END
