c=======================================================================
c
c     subroutine KERNELD                                      
c
c     Density Estimation by Non Parametric Gaussian Kernel method 
c
c-----------------------------------------------------------------------
      SUBROUTINE kerneld ( n, m, x, y, z, info )
c-----------------------------------------------------------------------
c
c     INPUT :
c            n      : number of histogram points (n>1)           integer
c            m      : number of grid points      (m>0)           integer
c            x      : returns/histogram (n)                       double
c            y      : grid (m)                                    double
c
c     OUTPUT :
c            z      : Gaussian kernel estimation (m)              double
c            info   : diagnostic argument                        integer
c
c-----------------------------------------------------------------------
c
      IMPLICIT NONE
c
c     arguments i/o
      INTEGER n, m, info
      DOUBLE PRECISION x(*), y(*), z(*)
c
c     local variables
      INTEGER i, j
      DOUBLE PRECISION h, hh, stdx, coef, pi, tmp, sumz
      PARAMETER ( pi = 3.1415926535898 )
c
c     external subroutines
      EXTERNAL VOLAT
c     
c     intrinsic functions
      INTRINSIC sqrt, exp
c
c-----------------------------------------------------------------------
c
c     initializations
      info = 0
      DO i = 1,m
        z(i) = 0.0
      ENDDO
c
c     standard deviation (volatility) of x
      CALL VOLAT ( n, x, stdx, info )
      IF (info .LT. 0) RETURN
c
c     optimal choice of h for Gaussian kernel (cf. D. Bosq and al. p85)
      h = (1.059*stdx)/(n**0.2)
c
c     coefficient 
      coef = 1./(SQRT(2.0*pi))
      hh   = 2.0*(h*h)
c
c     Gaussian kernel estimation at z
      sumz = 0.0
      DO i = 1,m
        tmp = 0.0
          DO j = 1,n
             tmp = tmp + coef * EXP((-(y(i)-x(j))**2)/(hh))
          ENDDO
          z(i) = tmp / (n*h)
      ENDDO
      RETURN
      END
